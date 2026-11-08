"""Exact consequences of the existing full-ground spin-projection theorem."""
from fractions import Fraction as F

from model.strong_coupling_refinement.bounds import finite_graph


def positive(x,name):
    if type(x)not in (int,F)or x<=0:raise ValueError(name+' must be a positive exact rational')
    return F(x)


def common_clock(coupling,contact,edges,spacing,clock,energy_budget=F(1)):
    c=positive(coupling,'hopping');u=positive(contact,'contact')
    p=positive(spacing,'spacing');q=positive(clock,'common clock')
    budget=positive(energy_budget,'energy budget')
    source=finite_graph(c,u,edges)
    if p>1:raise ValueError('spatial refinement requires p at most one')
    ratio=u/c;j=source['spin_exchange'];gap=source['charge_gap_lower']
    assert gap/j==ratio*(ratio/2-edges)
    assert gap/j>=ratio*ratio/4>=4*edges*edges
    rate=q*j*p
    return {'source_ratio_U_over_c':ratio,'spin_exchange':j,
            'unscaled_charge_gap_lower':gap,'ballistic_spin_clock_rate':rate,
            'rescaled_charge_gap_lower':q*gap,
            'charge_gap_per_ballistic_spin_rate':ratio*(ratio/2-edges)/p,
            'universal_gap_per_spin_rate_lower':4*F(edges*edges)/p,
            'complete_charged_weight_upper':min(F(1),budget/(q*gap)),
            'original_odd_creation_probability_upper_at_energy_budget':min(F(1),budget/(q*gap)),
            'spin_rate_upper_if_both_conjugate_gaps_at_most_budget':budget*p/(ratio*(ratio/2-edges)),
            'global_ground_projection_hypothesis_used':True,
            'local_dressed_spin_theorem_outside_domain_excluded':False}


def chemical_pair(plus_gap,minus_gap,shift,clock):
    plus=positive(plus_gap,'addition gap');minus=positive(minus_gap,'removal gap')
    q=positive(clock,'common clock')
    if type(shift)not in (int,F):raise ValueError('exact signed chemical shift required')
    nu=F(shift)
    return {'same_reference_addition_energy':q*(plus-nu),
            'same_reference_removal_energy':q*(minus+nu),
            'conjugate_sum':q*(plus+minus),
            'sum_independent_of_chemical_shift':True,
            'same_reference_stable_against_these_two_sectors':plus>=nu and minus>=-nu,
            'negative_excitation_reinterpreted_as_ground_gap':False}
