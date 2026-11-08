"""Selected-source regional energy tails and limits of energy-only state control."""
from fractions import Fraction as F
from functools import lru_cache
import sympy as s

from .algebra import spins
from .bounds import ground_energy_cube


def regional_second_moment(side):
    row = ground_energy_cube(side)
    volume,internal,crossing = row['sites'],row['internal_edges'],row['crossing_edges']
    first = F(crossing)
    second = F(min(internal*crossing,crossing*crossing))
    return {'side':side,'sites':volume,'crossing_bonds':crossing,
            'regional_excitation_first_moment_over_J_upper':first,
            'regional_excitation_second_moment_over_J2_upper':second,
            'selected_ground_equation_density_second_moment_over_J2_upper':F(crossing**2,volume**2),
            'best_density_second_moment_over_J2_upper':second/(volume**2),
            'requires_selected_state_limit_of_finite_source_grounds':True,
            'volume_limit_precedes_zero_selection_field':True,
            'full_regional_density_matrix_convergence_modulus':False}


def excitation_tail(side, threshold_over_J):
    threshold = F(threshold_over_J)
    if threshold <= 0:
        raise ValueError('positive regional excitation threshold required')
    row = regional_second_moment(side)
    first = row['regional_excitation_first_moment_over_J_upper']
    second = row['regional_excitation_second_moment_over_J2_upper']
    return min(F(1),first/threshold,second/(threshold*threshold))


def source_gap_control(side, field=F(1,4)):
    from model.native_ordered_ground.source import perron_bound
    if type(side) is not int or side < 2:
        raise ValueError('source cube side at least two required')
    field = F(field)
    sites = side**3
    edges = []
    index = lambda x,y,z:(x*side+y)*side+z
    for x in range(side):
        for y in range(side):
            for z in range(side):
                here = index(x,y,z)
                if x+1 < side: edges.append((here,index(x+1,y,z)))
                if y+1 < side: edges.append((here,index(x,y+1,z)))
                if z+1 < side: edges.append((here,index(x,y,z+1)))
    bound = perron_bound(sites,tuple(edges),F(1),field)
    gap = F(bound['Hamiltonian_gap_lower'])
    # Onsite selection terms change the regional spectral width, so use only
    # the first-moment boundary estimate, which is unchanged at finite field.
    tail = min(F(1), F(6*side*side)/gap)
    assert tail == 1
    return {'side':side,'selection_field_over_J':str(field),
            'actual_finite_Doeblin_gap_lower_over_J':str(gap),
            'regional_excited_subspace_probability_upper_from_boundary_and_gap':str(tail),
            'this_bound_certifies_a_pure_regional_ground':False,
            'Doeblin_bound_is_an_upper_bound_on_the_true_gap':False}


@lru_cache(maxsize=1)
def source_phase_control():
    from model.spin_charge_dispersion import dimer
    order = dimer.certificate()
    lower = F(order['simple_strict_magnetization_lower'])
    m = s.symbols('m',positive=True)
    sigma_x = s.Matrix([[0,1],[1,0]])
    plus = (s.eye(2)+2*m*sigma_x)/2
    minus = (s.eye(2)-2*m*sigma_x)/2
    difference = plus-minus
    assert difference*difference == 4*m*m*s.eye(2)
    spin = spins(2)
    h = sum((spin[0][a]*spin[1][a] for a in range(3)),s.zeros(4))-s.eye(4)/4
    gamma = 4*spin[0][1]*spin[1][1]
    assert gamma*h*gamma == h
    assert gamma*spin[0][0]*gamma == -spin[0][0]
    return {'states':'opposite maximal-order phases of the same zero-field source',
            'source_Hamiltonian_and_all_its_energy_statistics_equal':True,
            'one_site_half_trace_distance':'2m',
            'one_site_half_trace_distance_strict_lower':str(2*lower),
            'replayed_source_magnetization_bound':order,
            'both_states_are_the_same_plus_field_selected_state':False,
            'source_energy_convergence_implies_regional_state_convergence':False}


def regional_tomography_bound(pauli_differences):
    differences = tuple(F(x) for x in pauli_differences)
    # Returns the exact square of an upper bound on half trace distance.
    return min(F(1),sum((x*x for x in differences),F())/4)
