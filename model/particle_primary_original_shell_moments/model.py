"""Actual original-state off-diagonal first moments in the even density sector."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_energy_tightness import model as energy
from model.particle_primary_coherent_odd_response import kernel_reduction as kernels
from model.particle_primary_joint_spectroscopy import spectroscopy
from model.particle_primary_sector_threshold import model as sectors
from model.particle_primary_gravity_premises import model as modular
from model.particle_primary_occupation_bound import model as occupation
from model.particle_spatial_modular_response import regional

SOURCE_INPUTS=tuple(dict.fromkeys((
 'model/particle_primary_original_shell_moments/__init__.py',
 'model/particle_primary_original_shell_moments/__main__.py',
 'model/particle_primary_original_shell_moments/model.py',
 'model/particle_primary_original_shell_moments/test_model.py',
 'model/particle_primary_original_shell_moments/PROOF.md',
 'model/particle_primary_original_shell_moments/README.md',
 'model/particle_primary_coherent_odd_response/KERNEL_PROOF.md',
 'model/particle_spatial_modular_response/PROOF.md',
 'model/particle_primary_sector_threshold/PROOF.md',
 'model/particle_primary_gravity_premises/PROOF.md',
 'model/particle_primary_occupation_bound/PROOF.md',
 *energy.SOURCE_INPUTS,*spectroscopy.SOURCE_INPUTS,*regional.SOURCE_INPUTS,
)))


@lru_cache(maxsize=1)
def exact_geometry():
    cells=kernels.CELLS
    adjacency=s.Matrix([[int(kernels.distance(x,y)==1) for y in cells] for x in cells])
    W=s.Matrix([[(-1)**sum(a*b for a,b in zip(q,x)) for x in cells] for q in cells])/s.sqrt(8)
    laplacian=3*s.eye(8)-adjacency
    actual=energy.density_commutator_control()
    if W*laplacian*W.T!=s.diag(*(2*sum(q) for q in cells)):
        raise ArithmeticError('actual cube density shell Laplacian changed')
    orbit=kernels.orbit_control()
    hopping_symbol=s.Symbol('tau',positive=True)
    shells=(3*hopping_symbol/2,-hopping_symbol/2,0,0)
    if kernels.walsh_eigenresponses(shells)!=tuple(j*hopping_symbol for j in range(4)):
        raise ArithmeticError('original Ward shell moments differ from the spectral transform')
    if any(row['internal_dimension']<=0 for row in actual['complete_native_blocks']):
        raise ArithmeticError('a native source block was omitted')
    return {'source_Wilson_quantum_link_control':actual,'source_cube_pair_orbits':orbit,
        'ordered_pair_shell_sizes':[8,24,24,8],
        'density_first_moment_shell_coefficients_of_negative_bond_expectation':['3/2','-1/2','0','0'],
        'Walsh_eigenmoment_coefficients':[0,1,2,3],
        'Walsh_shell_multiplicities':[1,3,3,1],
        'source_matter_modes_per_cell':60,'centered_local_number_norm_upper':30,
        'complete_native_dimension_sum':sum(row['internal_dimension'] for row in actual['complete_native_blocks']),
        'full_physical_hopping_bond_norm_upper':30,
        'all_nonadjacent_off_diagonal_first_moments_exactly_zero':True}


def finite_record_error(norm,resolution,tail,implementation):
    """Original unnormalized signed-mean error for the complete two clocks."""
    norm,eta,tail,implementation=map(F,(norm,resolution,tail,implementation))
    if norm<=0 or eta<=0 or min(tail,implementation)<0:
        raise ValueError('positive source norm and resolution with nonnegative errors required')
    return 2*eta+56*norm*(2*tail+implementation)


@lru_cache(maxsize=1)
def native_sign_control():
    """Actual filled-beta anchor and all native first-action Haar/CAR norms."""
    rows=sectors.hopping_packets((0,0,0,0,0))
    anchor=sectors.local_singlet((0,0,0,0,0))
    if anchor['norm_squared']!=1 or anchor['defect_count']!=0:
        raise ArithmeticError('the filled-beta Gauss anchor changed')
    for row in rows:
        d=row['dimension']
        if row['unchanged_defect_packet_norm_squared'] or row['two_added_defect_packet_norm_squared']!=12*d:
            raise ArithmeticError('native first-action physical anchor norms changed')
        if any(axis[2]!=d or axis[0] for axis in row['affected_axis_norms']):
            raise ArithmeticError('one of the actual Wilson-axis packet norms changed')
    cells=kernels.CELLS
    bonds=[(x,y) for x in cells for y in cells if x<y and kernels.distance(x,y)==1]
    if len(bonds)!=12 or any((-1)**(sum(x)+sum(y))!=-1 for x,y in bonds):
        raise ArithmeticError('native sublattice parity does not reverse every actual bond')
    return {'physical_filled_beta_anchor':{k:v for k,v in anchor.items() if k!='state'},
        'all_native_first_action_packets':rows,'actual_bipartite_bonds_checked':len(bonds),
        'per_native_anchor_variance_factor_of_dimension':12,
        'per_native_physical_Hilbert_Schmidt_norm_squared_lower_factor_of_dimension':24,
        'anchor_used_as_original_input_state':False}


def native_sign_bound(member,floor):
    """Strict native hopping signs in the actual selected polynomial state."""
    control=native_sign_control();slope=modular.selected_modular_slope(member)
    beta=F(slope['modular_slope_lower']);V=member['source']['side']**3
    if beta<=0:raise ArithmeticError('strict selected polynomial slope is required')
    rows=[]
    for row in control['all_native_first_action_packets']:
        d=row['dimension']
        factored={**floor,'operator':'positive scalar bound'}
        rows.append({'Q6':row['Q6'],'native_dimension':d,
            'density':'n_(r,x)=N_(r,x)-2d I','density_norm_upper':2*d,
            'physical_hopping_HS_norm_squared_lower':24*d,
            'negative_native_bond_expectation_lower_factored':{**factored,'rational_multiplier':2*d*beta},
            'four_original_first_moment_factors_of_tau_r':[F(j,4*d*d*V) for j in range(4)],
            'nonconstant_smallest_first_moment_lower_factored':{**factored,'rational_multiplier':beta/(2*d*V)},
            'native_density_first_moment_shells_strictly_nonconstant':True})
    return {'original_selected_polynomial_slope':slope,'source_physical_density_floor_factored':floor,
        'actual_native_control':control,'native_density_channels':rows,
        'exact_native_identity':'U_r H U_r=H-2T_r, U_r=(-1)^(sum_black N_(r,x))',
        'every_native_total_hopping_expectation_strictly_negative':True,
        'finite_total_density_clock_resolves_each_native_lower_bound':False,
        'native_density_channels_identified_with_particle_species':False}


def density_success_bound(member,tau_lower,mean_error):
    """Complete event-weight cap and separately labelled conditional statistic."""
    inherited=occupation.bind_member(member);occupation.verify_bound(member,inherited)
    D=F(inherited['result']['defect_count_upper']);V=member['source']['side']**3
    success=min(F(1),D/(30*V))
    original=tau_lower/600;reported=original-7*mean_error
    if success<=0 or reported<=0:
        raise ArithmeticError('source density-event mean or success bound failed')
    return {'occupation_receipt_sha256':contract.digest(inherited),
        'source_defect_count_upper':D,'local_full_Fock_inequality':'(N_x-30I)^2 <=30 D_x',
        'original_complete_density_success_weight_upper':success,
        'same_success_weight_with_complete_source_clocks':True,
        'all_eight_original_signed_means_sum_lower':original,
        'all_eight_reported_signed_means_sum_lower':reported,
        'conditional_on_density_success_ideal_signed_mean_lower':original/success,
        'conditional_on_density_success_reported_signed_mean_lower':reported/success,
        'conditional_statistic_is_original_event_weight':False,
        'conditional_statistic_is_isolated_particle_energy_or_mass':False}


@lru_cache(maxsize=1)
def exact_spectral_countercheck():
    """Independent exact finite spectral sum checks the Ward factors and signs."""
    cells=kernels.CELLS;V=8
    A=s.Matrix([[int(kernels.distance(x,y)==1) for y in cells] for x in cells]);H=-A
    W=s.Matrix([[(-1)**sum(a*b for a,b in zip(q,x)) for x in cells] for q in cells])/s.sqrt(8)
    energies=[2*sum(q)-3 for q in cells];Q=4*s.eye(V)-H;rho=Q*Q/s.trace(Q*Q)
    diagonal=W*rho*W.T
    spectral=s.zeros(V);fields=[]
    for x in range(V):
        field=s.zeros(V);field[x,x]=1;fields.append(W*field*W.T)
    for x,y in product(range(V),repeat=2):
        spectral[x,y]=sum(diagonal[i,i]*(energies[f]-energies[i])*fields[x][f,i]*fields[y][f,i]
                          for i,f in product(range(V),repeat=2))
    t=s.trace(rho*(-s.eye(V)[:,0]*s.eye(V)[:,1].T-s.eye(V)[:,1]*s.eye(V)[:,0].T))
    expected=(-t/2)*(3*s.eye(V)-A)
    if spectral!=expected:raise ArithmeticError('direct original spectral sums differ from the Ward shell matrix')
    diagonalized=W*spectral*W.T
    if diagonalized!=s.diag(*(-t*sum(q) for q in cells)):
        raise ArithmeticError('independent energy eigenbasis lost shell ratios')
    return {'domain':'eight-state one-particle cube control, not the primary source',
        'density':'(4I-H)^2/Tr((4I-H)^2)','negative_bond_expectation':F(-t),
        'four_first_moments':[F(-t*j) for j in range(4)],
        'all64_original_spectral_kernel_entries_checked':True,
        'control_preparation_replaces_primary_polynomial_density':False}


def bind_member(member):
    contract.verify_member(member)
    if member['source']['level']!=0:raise ValueError('the four-shell original density calculation uses member zero')
    controls=exact_geometry();inherited=energy.selected_density_sum_rule(member)
    V=member['source']['side']**3;B=member['source']['interior_hopping_bonds']
    if V!=8 or B!=12 or F(member['source']['spacing'])!=1:
        raise ArithmeticError('actual member-zero graph changed')
    tau_lower=F(inherited['selected_negative_hopping_expectation_lower'])/B
    tau_upper=F(controls['full_physical_hopping_bond_norm_upper'])
    first_unit=tau_lower/(900*V)
    if not 0<tau_lower<=tau_upper or first_unit!=F(inherited['original_signed_first_moment_lower']):
        raise ArithmeticError('original first-shell moment fails the authenticated source sum rule')
    source=member['preparation_plan']['source']
    norm=max(abs(F(source['spectral_lower'])),abs(F(source['spectral_upper'])))
    eta=first_unit/16;budget=first_unit/(1024*norm)
    clock=spectroscopy.clock_plan(member,eta,budget,budget)
    tail=F(clock['clock_plan']['complete_clock_tail_upper']);implementation=F(clock['two_complete_clock_error_upper'])
    error=finite_record_error(norm,eta,tail,implementation)
    if not error<first_unit/2 or not first_unit-2*error>0:
        raise ArithmeticError('same-source complete clocks do not resolve the first-moment contrast')
    shell_rows=[]
    for weight,multiplicity in enumerate(kernels.MULTIPLICITIES):
        ideal=(F(weight)*first_unit,F(weight)*tau_upper/(900*V))
        reported=(F(),F()) if weight==0 else (ideal[0]-error,ideal[1]+error)
        shell_rows.append({'Walsh_weight':weight,'multiplicity':multiplicity,
            'exact_first_moment_factor_of_tau':F(weight,900*V),
            'original_signed_transfer_first_moment_interval':ideal,
            'reported_original_signed_transfer_first_moment_interval':reported,
            'elastic_ideal_measure_only':weight==0,
            'positive_original_inelastic_energy_weight_proved':weight>0})
    floor=regional.common_contract_payload(member,1)['physical_density_floor_factored']
    result=contract.serial({'member_id':member['member_id'],**member['ids'],
        'source_input_is_original_selected_density':True,
        'probe':'n_x=N_x-30I, using every one of the60 native CAR modes',
        'probe_CAR_parity':0,'probe_native_occupation_increment':[0,0,0,0,0],
        'odd_G1_H_trion_shell_values_computed':False,
        'same_source_density_sum_rule_sha256':contract.digest(contract.serial(inherited)),
        'source_controls':controls,'negative_bond_expectation_interval':[tau_lower,tau_upper],
        'first_moment_shell_matrix':'M^(1)=tau*(3I-Adj_cube)/2, tau=-Tr(rho H_hop)/12>0',
        'first_moment_signed_shell_intervals':[[3*tau_lower/2,3*tau_upper/2],[-tau_upper/2,-tau_lower/2],[0,0],[0,0]],
        'full_menu_success_Kraus':'sum_x (-1)^(q dot x) n_x/(30*8)',
        'full_menu_failure_Kraus':'sqrt(I-(n_x/30)^2)/sqrt(8), with every address retained',
        'spatial_success_records':8,'address_failure_records':8,'records_before_source_clocks':16,
        'per_source_field_normalization':900*V,'all_failures_and_reference_retained':True,
        'spatial_shell_records':shell_rows,'complete_two_energy_clock_plan':clock,
        'finite_signed_mean_error_upper':error,
        'all_nonzero_shell_reported_signed_means_strictly_positive':True,
        'neighboring_nonzero_shell_reported_mean_difference_lower':first_unit-2*error,
        'zero_shell_reported_signed_mean_exactly_zero':True,
        'zero_shell_source_success_weight_lower_factored':floor,
        'native_hopping_sign_and_density_shells':native_sign_bound(member,floor),
        'complete_density_success_and_conditional_statistics':density_success_bound(member,tau_lower,error),
        'zero_shell_positive_weight_witness':'empty Fock and constant-link vector gives K_0=-I on that vector',
        'zero_and_nonzero_shell_normalized_ideal_energy_shapes_differ':True,
        'nonzero_shell_normalized_energy_shapes_proved_distinct':False,
        'complete_clock_histogram_numerically_executed':False,
        'conserved_momentum_or_relativistic_dispersion_proved':False,
        'isolated_particle_branch_or_uniform_refinement_residue_proved':False,
        'physical_mass_table_GeV':[]})
    return {'schema':'oph.primary_original_shell_moments.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,receipt):
    if receipt!=bind_member(member):raise ValueError('original spectral shells, state or clock changed')


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_original_shell_moments_certificate.v1',
        'independent_spectral_sum_control':contract.serial(exact_spectral_countercheck()),
        'selected_member_bound':bind_member(contract.member(0)),'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('source changed during original shell moment replay')
    return result
