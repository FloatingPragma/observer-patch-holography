"""A dynamical signed extension and exact complete-input preparation work."""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import comb
import numpy as np
from scipy.linalg import eigvalsh
from model.openwave_nuclear_em_transfer.source import bound,A,B
from model.openwave_light_nuclear_curve.source import PARAMETERS,UNIT
from model.openwave_light_nuclear_curve.certificate import INPUTS,singlet_dimension
from model.openwave_capacity_stress.source import enclosure,sqrt_interval
from model.coherent_defect_wave.chiral_joint.model import fermion

OMEGA=A
COUPLING=PARAMETERS['t']

def shift(coupling=COUPLING,omega=OMEGA):
    coupling,omega=F(coupling),F(omega)
    if omega<=0:raise ValueError('strictly positive oscillator stiffness required')
    return coupling*coupling/(2*omega)

def four_site_occupancy(mask,cells):
    return sum(all(((mask>>(12*x))&4095).bit_count()==3 for x in range(start,start+4))
               for start in range(max(0,cells-3)))

def noncommuting_native_control():
    # Each local triple contains all three colours of one internal orbital.
    mask=sum(1<<(12*x+j)for x in range(4)for j in((0,1,2)if x%2==0 else(3,4,5)))
    a=fermion(mask,0,False);assert a is not None
    b=fermion(a[0],12,True);assert b is not None
    assert four_site_occupancy(mask,4)==1 and four_site_occupancy(b[0],4)==0
    return {'native_cells':4,'fixed_quark_number':12,'initial_occupations':[3,3,3,3],
            'one_native_hopping_successor_occupations':[2,4,3,3],
            'nonzero_O_H_hopping_commutator_coefficient_exact':str(PARAMETERS['t']*a[1]*b[1]),
            'static_elimination_exact_for_full_dynamics':False}

def casimir_charge_second(cells,suffix):
    L,r=cells,suffix
    if L<1 or not 0<=r<=L:raise ValueError('invalid complete singlet carrier or suffix')
    gamma=F(9*L*(4*L+3),4*(16*L*L-1))
    return F(r*r,4)+gamma*(F(10*r,9)-F(r*r,9*L))

def full_uniform_work(cells):
    gamma=F(9*cells*(4*cells+3),4*(16*cells*cells-1))
    suffix=[casimir_charge_second(cells,r)for r in range(1,cells)]
    electric=9*A*sum(suffix,F(0))/2
    return {'cells':cells,'complete_fixed_N_colour_singlet_dimension':singlet_dimension(cells),
            'input':'Normalized identity on the entire N=3L global colour-singlet carrier, with every charge block and compensating reference.',
            'orbital_irrep':'U(4L) Young diagram (3^L)',
            'quadratic_trace_coefficient_exact':str(gamma),
            'all_suffix_charge_square_means_exact':list(map(str,suffix)),
            'bare_Gauss_electric_mean_source_exact':str(electric),
            'full_field_service_mean_work_source_interval':[str(electric),str(electric+B)],
            'neutral_oscillator_vacuum_initial_service_work_exact':'0',
            'conditional_fixed_charge_quartet_mean_substituted':False}

def character_control():
    # Direct Schur-character monomials, including every weight multiplicity.
    rows=list(combinations_with_replacement(range(8),3))
    charge=[F(2,3),F(2,3),F(-1,3),F(-1,3),F(0),F(0),F(0),F(0)]
    values=[sum((charge[k]for k in left+right),F(0))for left in rows for right in rows
            if all(x<y for x,y in zip(left,right))]
    assert len(values)==singlet_dimension(2)==2520
    mean=sum(values,F(0))/len(values);second=sum((q*q for q in values),F(0))/len(values)
    assert mean==F(1,2) and second==casimir_charge_second(2,1)==F(68,63)
    return {'method':'Complete semistandard-tableau character of shape (3,3) on eight orbitals, independently enumerated.',
            'complete_weight_monomials':len(values),'native_suffix_charge_mean_exact':str(mean),
            'native_suffix_charge_square_mean_exact':str(second),'all_weight_multiplicities_retained':True}

def oscillator_controls():
    rows=[]
    for cutoff in(16,32):
        annihilator=np.diag(np.sqrt(np.arange(1,cutoff+1,dtype=float)),1)
        H=float(OMEGA)*np.diag(np.arange(cutoff+1,dtype=float))
        H+=float(COUPLING)/np.sqrt(2)*(annihilator+annihilator.T)
        eigenvalue=float(eigvalsh(H,subset_by_index=[0,0])[0])
        assert abs(eigenvalue+float(shift()))<1e-6
        rows.append({'number_cutoff_for_control_only':cutoff,'lowest_eigenvalue_diagnostic':eigenvalue,
                     'distance_to_exact_full_oscillator_ground_diagnostic':abs(eigenvalue+float(shift()))})
    assert rows[1]['distance_to_exact_full_oscillator_ground_diagnostic']<2e-14
    return rows

def nuclear(cells,charge,name,coupling=COUPLING):
    parent=bound(cells,charge);blocks=max(0,cells-3);decrease=blocks*shift(coupling)
    original=tuple(map(F,parent['full_dynamic_ground_GeV_interval']))
    result=tuple(x-decrease*UNIT for x in original)
    target=F(INPUTS['nuclear_rest_energy_GeV'][name])
    upper_rows=[]
    for row in parent['all_bare_flux_zero_spin_multiplets']:
        upper_rows.append({'spin':row['spin'],'rank':row['rank'],
            'all_displaced_quartet_variational_levels_source_intervals':[
                [str(F(lo)+B-decrease),str(F(hi)+B-decrease)]
                for lo,hi in row['bare_cycle_flux_zero_variational_eigenvalue_source_intervals']]})
    return {'name':name,'cells':cells,'charge':charge,'local_four_site_field_modes':blocks,
        'full_ground_GeV_interval':list(map(str,result)),'full_ground_GeV_display':enclosure(result),
        'parent_full_ground_GeV_interval':parent['full_dynamic_ground_GeV_interval'],
        'maximum_field_induced_ground_lowering_GeV_exact':str(decrease*UNIT),
        'target_GeV':str(target),'target_role':'inherited calibration'if name=='triton'else'held-out comparison',
        'target_in_certified_interval':result[0]<=target<=result[1],
        'comparison':'excluded_below'if target<result[0]else'excluded_above'if target>result[1]else'not_resolved_by_this_enclosure',
        'all_quartet_variational_spin_blocks':upper_rows,
        'variational_levels_are_complete_source_excited_eigenvalues':False,
        'full_matter_EM_flux_and_oscillator_complements_retained':True,
        'charge_quartet_preparation_success_exact':parent['quartet_filter_success_on_parent_uniform_input_exact'],
        'charge_quartet_preparation_failure_exact':parent['quartet_filter_failure_exact'],
        'conditional_oscillator_free_energy_after_displacement_source_exact':str(decrease),
        'conditional_coupling_energy_after_displacement_source_exact':str(-2*decrease),
        'conditional_displacement_work_for_full_coupled_H_source_exact':str(-decrease),
        'nonselective_displacement_work_on_uniform_input_source_exact':str(-decrease*F(parent['quartet_filter_success_on_parent_uniform_input_exact']))}

def run():
    assert OMEGA==F(1,16) and COUPLING==F(1,8) and shift()==F(1,8)
    controls=character_control();commutator=noncommuting_native_control();oscillators=oscillator_controls()
    rows=[nuclear(L,Z,name)for name,L,Z in(('triton',3,1),('helion',3,2),('alpha',4,2))]
    assert rows[0]['full_ground_GeV_interval']==rows[0]['parent_full_ground_GeV_interval']
    assert rows[1]['comparison']=='excluded_above' and rows[2]['comparison']=='excluded_below'
    alpha=rows[2];gap=F(alpha['parent_full_ground_GeV_interval'][0])-F(alpha['target_GeV'])
    necessary_square=2*OMEGA*gap/UNIT;needed=sqrt_interval(necessary_square)
    assert needed[0]>COUPLING
    return {'schema':'oph.native_nuclear_dynamical_four_body_attraction.v1',
        'verdict':'PASS_STABLE_FULL_FIELD_ATTRACTION_COMPLETE_UNIFORM_WORK_AND_HELD_OUT_EXCLUSIONS',
        'declared_action':{'parent':'Unchanged native nuclear strong source and full Gauss-cycle electromagnetic extension.',
            'local_operator':'O_b=product over four consecutive distinct sites of 1_(N_x=3); each O_b is a bounded gauge-invariant projector.',
            'neutral_field':'One complete L2(R) oscillator per four-site block, with [q_b,p_b]=i.',
            'Hamiltonian':'H_EM + sum_b omega/2*(p_b^2+q_b^2-1) + g*q_b*O_b.',
            'omega_source_exact':str(OMEGA),'g_source_exact':str(COUPLING),
            'coefficient_ancestry':'omega is the inherited electric architecture coefficient1/16; g is the inherited native hopping coefficient1/8. These are declared scalar-coupling choices, with no nuclear target used to choose them.',
            'zero_on_all_N_less_than_12':True,'electron_light_calibration_and_triton_coefficients_changed':False,
            'same_existing_source_without_extension_claimed':False},
        'complete_square':{'identity':'omega*a_b^*a_b + g*q_b O_b = omega*(a_b^*+g O_b/(sqrt(2)omega))*(a_b+g O_b/(sqrt(2)omega)) - g^2 O_b^2/(2omega).',
            'full_operator_lower_bound':'H >= H_EM - [g^2/(2omega)] sum_b O_b >= H_EM - [g^2/(2omega)] number_of_blocks.',
            'stability_bound_extensive_in_blocks':True,'signed_coupling_escapes_positive_extension_hypothesis':True,
            'dynamical_oscillator_replaced_by_a_static_potential':False,'native_noncommutation_control':commutator,
            'full_scalar_ground_for_fixed_O1_source_exact':str(-shift()),'finite_Fock_diagnostics_only':oscillators},
        'complete_ground_bounds':rows,
        'full_uniform_preparation_work':[full_uniform_work(L)for L in(3,4)],'independent_character_control':controls,
        'necessary_alpha_strength':{'necessary_g_squared_source_lower_exact':str(necessary_square),
            'necessary_absolute_g_source_interval':list(map(str,needed)),
            'necessary_absolute_g_source_display':enclosure(needed),
            'declared_g_source_exact':str(COUPLING),'bound_used_to_choose_source_coefficient':False,
            'sufficient_to_match_alpha_or_helion':False},
        'complete_observer_records':{'initial':'Full fixed-N singlet identity, Gauss-dressed cycle n=0, all charge/reference blocks, and neutral oscillator vacua.',
            'filter':'Complete fixed-charge native quartet projector and its full failure; zero displacement on failure.',
            'success_field_trial':'On every quartet-success state O_b=1, apply the normalized full coherent displacement -g/omega; its infinite oscillator tail is exact.',
            'displacement_work':'The free oscillator receives+c per block and its coupling energy changes by-2c, so full-source work is-c on the success branch; failure has zero displacement work.',
            'terminal':'Full coupled-H energy spectral measure with every energy complement; nonselective terminal energy work0 and arbitrary idle reference.',
            'finite_local_compilation_of_exact_spectral_service':False},
        'held_out_helion_or_alpha_reproduced':False,'complete_breakup_or_halo_thresholds_proved':False,
        'new_empirical_fit':False,'new_match_with_nature':False}
