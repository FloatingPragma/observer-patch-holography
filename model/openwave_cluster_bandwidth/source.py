"""Cluster-additive native strange bands with an open two-particle energy shell."""
from fractions import Fraction as F
from pathlib import Path
import json
from model.openwave_hadron_refinement_flow.positivity import joint_stability
from model.openwave_baryon_dispersion.source import phase_record,repair_matrix
from model.openwave_capacity_stress.source import enclosure
from .lambda_source import ground as lambda_ground
from .meson import ground as meson_ground, native_two_cell_check, data as meson_data, degree_ancestry


def minimum_pair(nq,na):
    return min(nq,na)


def local_completion(counts):
    sites=sum(minimum_pair(*row)for row in counts)
    bonds=sum(minimum_pair(a[0]+b[0],a[1]+b[1])-minimum_pair(*a)-minimum_pair(*b)
              for a,b in zip(counts,counts[1:]))
    return sites+bonds


def angle(t,terms=32):
    t=F(t)
    if not 0<=t<=F(1,2):raise ValueError('arctangent certification domain')
    value=sum((-1)**n*t**(2*n+1)/F(2*n+1)for n in range(terms))
    next_term=(-1)**terms*t**(2*terms+1)/F(2*terms+1)
    return min(2*value,2*(value+next_term)),max(2*value,2*(value+next_term))


def shell(t):
    t=F(t);c=(1-t*t)/(1+t*t);s=2*t/(1+t*t)
    pro=phase_record(c,s,True,True);pion=meson_ground(c,-s)
    a=tuple(map(F,pro['all_five_full_source_GeV_intervals'][0]));b=tuple(map(F,pion['full_ground_GeV_interval']))
    pair=(a[0]+b[0],a[1]+b[1])
    return {'half_angle_tangent':str(t),'phase_interval':enclosure(angle(t)),
       'proton':pro,'pion':pion,'complete_pair_energy_GeV_interval':[str(x)for x in pair],
       'pair_energy_GeV_display':enclosure(pair)}


def run():
    lam=lambda_ground();a,b=map(F,lam['full_Lambda_ground_GeV_interval'])
    left,right=shell(F(1,6)),shell(F(1,5))
    lower_pair=tuple(map(F,left['complete_pair_energy_GeV_interval']));upper_pair=tuple(map(F,right['complete_pair_energy_GeV_interval']))
    assert lower_pair[1]<a and upper_pair[0]>b
    base=joint_stability();c=base['C_cost_per_pair'];floor=base['uniform_floor']-c/2
    assert floor>0
    controls=[]
    for counts in (((1,0),(0,1),(1,0)),((2,0),(0,0),(0,1)),((2,0),(0,0),(0,0),(0,1))):
        value=local_completion(counts);N=sum(sum(row)for row in counts);assert 0<=value<=N
        controls.append({'counts':[list(row) for row in counts],'finite_range_completion':value,'global_minimum_control':min(sum(q for q,a in counts),sum(a for q,a in counts))})
    from model.openwave_baryon_dispersion.source import ROOT
    units=json.loads((ROOT/'model/openwave_string_positive_cone/results.json').read_text())['common_length_attachment']
    hc=F(json.loads((ROOT/'model/openwave_transverse_gravity/inputs.json').read_text())['hbar_c_GeV_fm'])
    scale=hc/F(units['a_fm']);momentum=(angle(F(1,6))[0]*scale,angle(F(1,5))[1]*scale)
    zero_pion=meson_ground(F(1),F(0));zero_kaon=meson_ground(F(1),F(0),True)
    from model.openwave_hadron_refinement_flow.certified_rest import sqrt_upper
    from model.openwave_hadron_refinement_flow.analytic import constants
    _,_,_,p,_,unit=constants(*__import__('model.openwave_hadron_refinement_flow.parameters',fromlist=['degree_coefficients']).degree_coefficients())
    oldp=phase_record(F(1),F(0),False,True)
    old_upper=(F(oldp['all_five_compact_source_intervals'][0][1])+2*sqrt_upper(3)*p['t']+meson_data(False)[0])*unit
    assert old_upper<a
    return {'schema':'oph.cluster_additive_strange_bandwidth.v1',
       'verdict':'PASS_POSITIVE_LOCAL_STRONG_SOURCE_AND_CERTIFIED_OPEN_P_PI_ENERGY_SHELL',
       'finite_range_source':{'construction':'Möbius extension of the complete native two-cell strong operator; one-site vacuum restrictions plus connected nearest-neighbour terms and the declared Gauss electric energy.',
          'complete_two_cell_operator_preserved':True,'strange_mass_parameters_refitted':False,
          'C_completion':'sum_x min(nq_x,na_x)+sum_edges[min(nq_x+nq_y,na_x+na_y)-min(nq_x,na_x)-min(nq_y,na_y)]',
          'C_completion_bound':'0<=M<=total_number on a path',
          'global_minimum_is_used':False,'separated_cluster_interactions_vanish':True,
          'full_Fock_number_floor_source_exact':str(floor),'count_controls':controls,
          'correlated_propagation':'One positive inverse-native-rest-matrix rule for every triquark and quark-antiquark compact channel; full Fock complements retained.',
          'one_species_selected_by_projector':False},
       'Lambda_rest':lam,'pion_rest':zero_pion,'kaon_rest':zero_kaon,
       'native_meson_intertwiners':[native_two_cell_check(False),native_two_cell_check(True)],
       'meson_degree_ancestry':degree_ancestry(),
       'energy_shell_bracket':{'below':left,'above':right,
          'below_margin_GeV_exact':str(a-lower_pair[1]),'above_margin_GeV_exact':str(upper_pair[0]-b),
          'opposite_momenta_GeV_interval':enclosure(momentum),
          'argument':'Continuous isolated full ground bands and strict enclosing inequalities imply an energy-conserving p+pi shell for some momentum in the bracket.'},
       'locality_only_control':{'all_momenta_p_pi_ground_energy_upper_GeV':str(old_upper),
          'Lambda_rest_lower_GeV':str(a),'insufficient_bandwidth':True,
          'bound':'Baryon compact phase-norm bound plus the complete meson co-located variational energy; no correlated kinetic term.'},
       'native_weak_scope':{'parent_nonzero_native_vertex_retained_as_control':True,
          'parent_two_cell_closed_channel_retained_as_control':True,
          'new_on_shell_weak_matrix_element_certified':False,
          'complete_outgoing_weak_decay_or_lifetime_proved':False,
          'missing_transfer':'The unfiltered native current requires a common scattering-state and interacting-vacuum transfer. Threshold opening alone does not execute a decay.'},
       'all_parent_two_cell_mass_fits_are_bulk_poles':False,
       'bulk_energies_role':'Consequences of inherited fitted coefficients and a declared finite-range/kinetic extension; the small bulk rest shifts are retained.'}
