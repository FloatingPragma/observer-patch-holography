"""Exact Gauss electric operator, complete variational bands and no-lowering bound."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import numpy as np
import sympy as sy
from model.openwave_light_nuclear_curve.source import projected_source,PARAMETERS,UNIT,lift
from model.openwave_light_nuclear_curve.certificate import INPUTS,graph
from model.openwave_light_nuclear_curve.exact import rational_matrix
from model.openwave_hadron_refinement_flow.certified_rest import eigen_intervals
from model.openwave_capacity_stress.source import enclosure

from model.openwave_transverse_completion.source import ELECTRIC as A, MAGNETIC as B


def gauss_flux(charges,n=0):
    charges=tuple(F(q)for q in charges);assert all((3*q).denominator==1 for q in charges)
    flux=[F(n)-3*sum(charges[e+1:])for e in range(len(charges)-1)]+[F(n)]
    reference=-3*sum(charges)
    assert all(flux[e]-flux[e-1]==3*charges[e]+(reference if e==0 else 0)for e in range(len(charges)))
    return flux,reference


def electric(charges,n=0):
    flux,_=gauss_flux(charges,n)
    return A*sum(x*x for x in flux)/2


@lru_cache(None)
def bound(cells,charge):
    d=projected_source(cells,F(INPUTS['three_body_coefficient']),True)
    parent=graph(cells,F(INPUTS['three_body_coefficient']));Iz=F(charge)-F(cells,2)
    refIz=F(-1,2)if cells%2 else F(0);shift=-PARAMETERS['delta']*(Iz-refIz)
    # Every low parent multiplet and its full complement enter the lower floor.
    inherited=min(F(lo)for row in parent['low_multiplets']for lo,hi in row['full_source_GeV_intervals'])+shift*UNIT
    Q=[np.diag(lift(np.diag([1,0,1,0]),x,cells)).real for x in range(cells)]
    potentials=[electric([int(Q[x][j])for x in range(cells)])for j in range(d['dimension'])]
    blocks=[]
    for twice_s in range(cells%2,cells+1,2):
        spin=F(twice_s,2)
        ids=np.flatnonzero((np.diag(d['Sz'])==float(spin))&(np.diag(d['Iz'])==float(Iz)))
        S=rational_matrix(d['S2'][ids][:,ids]);columns=(S-sy.Rational(spin*(spin+1))*sy.eye(len(ids))).nullspace()
        if not columns:continue
        P=sy.Matrix.hstack(*columns);G=P.T*P;H=sy.zeros(len(ids))
        for coefficient,matrix in d['rational_terms']:
            H+=sy.Rational(coefficient)*rational_matrix(matrix[ids][:,ids])
        H+=sy.diag(*[sy.Rational(potentials[j])for j in ids])
        intervals=eigen_intervals(P.T*H*P,G)
        blocks.append({'spin':str(spin),'rank':len(columns),
            'bare_cycle_flux_zero_variational_eigenvalue_source_intervals':[[str(lo),str(hi)]for lo,hi in intervals]})
    low=min(F(row['bare_cycle_flux_zero_variational_eigenvalue_source_intervals'][0][1])for row in blocks)
    upper=(low+B)*UNIT # <1-cos(loop)>_n=0 = 1; subtract E_gamma>=0.
    assert inherited<upper
    active=[potentials[j]for j in range(d['dimension'])if np.diag(d['Iz'])[j]==float(Iz)]
    electric_mean=sum(active,F(0))/len(active)
    mean_formula=9*A/2*(F(charge*(cells-1),2)+F(charge*(charge-1)*(cells-2),3))
    assert electric_mean==mean_formula
    rank=__import__('math').comb(cells,charge)*2**cells
    assert sum(row['rank']*(2*F(row['spin'])+1)for row in blocks)==rank
    success=F(rank,parent['complete_physical_N3A_colour_singlet_dimension'])
    return {'cells':cells,'charge':charge,'full_dynamic_ground_GeV_interval':[str(inherited),str(upper)],
      'full_dynamic_ground_GeV_display':enclosure((inherited,upper)),
      'bare_flux_zero_operator_norm_on_charge_quartet_source_exact':str(max(active)),
      'bare_flux_zero_operator_norm_on_charge_quartet_GeV_exact':str(max(active)*UNIT),
      'uniform_charge_quartet_electric_mean_source_exact':str(electric_mean),
      'uniform_charge_quartet_field_preparation_energy_above_vacuum_source_interval':[str(electric_mean),str(electric_mean+B)],
      'all_bare_flux_zero_spin_multiplets':blocks,'complete_bare_flux_zero_magnetic_dimension':rank,
      'complete_parent_high_Fock_complement_retained':True,'all_integer_cycle_fluxes_retained':True,
      'quartet_filter_success_on_parent_uniform_input_exact':str(success),
      'quartet_filter_failure_exact':str(1-success),
      'entire_cycle_flux_support_for_source':'all n in Z; n=0 specifies only the variational trial and complete preparation record.'}


def run():
    # Charge constraints are checked also on fractional and negative occupations.
    checks=0
    for charges in product((F(-2,3),F(-1,3),F(0),F(1,3),F(2,3)),repeat=3):
        for n in(-2,0,3):gauss_flux(charges,n);checks+=1
    rows=[]
    for name,L,Z in(('triton',3,1),('helion',3,2),('alpha',4,2)):
        row=bound(L,Z);target=F(INPUTS['nuclear_rest_energy_GeV'][name]);lo,hi=map(F,row['full_dynamic_ground_GeV_interval'])
        row.update(name=name,target_GeV=str(target),target_role='inherited calibration'if name=='triton'else'held-out comparison',
                   target_in_certified_interval=lo<=target<=hi,
                   comparison='excluded_below'if target<lo else'excluded_above'if target>hi else'not_resolved_by_this_enclosure')
        rows.append(row)
    assert rows[2]['comparison']=='excluded_below'
    alpha_gap=F(rows[2]['full_dynamic_ground_GeV_interval'][0])-F(rows[2]['target_GeV'])
    return {'schema':'oph.native_nuclear_compact_electromagnetic_transfer.v1',
      'verdict':'PASS_EXACT_DYNAMIC_GAUSS_EXTENSION_AND_HELD_OUT_ALPHA_OBSTRUCTION',
      'declared_extension':{'native_nuclear_source':'Same complete fixed-N colour-singlet three/four-cell source, all rational strong and triton coefficients unchanged.',
        'spatial_graph':'Cycle on the same cells; the original open-chain native terms are dressed along its spanning tree.',
        'field_H':'(a/2) sum E_edge^2 + b(1-Re U_loop), a=1/16,b=1/8, with full integer rotors.',
        'field_coefficients_from':'Exact architecture coefficients of openwave_transverse_completion; no alpha or nuclear target selects them.',
        'full_operator':'W (H_native tensor I_cycle) W* + H_field - E_gamma.',
        'W':'Exact Gauss isometry on every native Fock occupation and every integer cycle flux, with R=-3Q_total.',
        'external_compensating_reference_Hamiltonian':'0, a declared apparatus term; its charge record is retained.',
        'additional_external_positive_hopping':'Not added; all native nuclear hopping and charge-exchange terms receive the exact same tree dressing.',
        'true_field_vacuum_subtracted':True,'vacuum_energy_source_interval':['0',str(B)],
        'same_existing_192_mode_electromagnetic_source_claimed':False,
        'electron_unit_GeV_exact':str(UNIT),'electron_anchor_recalibrated':False,
        'all_parent_mass_equalities_preserved_by_extension_proved':False,
        'empirical_alpha_is_source_coefficient':False},
      'exact_Gauss_checks':checks,
      'effective_operator':{'zero_cycle_flux_electric_operator':'(9a/2) sum_(e=0)^(L-2) (sum_(x>e) Q_x)^2, Q_x=N_x/6+I_z,x.',
        'cycle_flux_generalization':'(a/2) sum_e [n-3 sum_(x>e) Q_x]^2 with final edge E=n.',
        'zero_flux_compression_magnetic_mean_exact':str(B),
        'zero_flux_to_full_flux_coupling_norm_exact':str(B/sy.sqrt(2)),
        'full_retardation_and_polarization':'The magnetic shifts n to n+/-1 are retained in the source; no static elimination is used for the full ground lower bound.',
        'substituted_1_over_distance_Coulomb_kernel':False,
        'P_to_high_matter_Q_electric_coupling_exact':'0; native charge is diagonal and preserves every local colour/spin block.'},
      'complete_ground_bounds':rows,
      'positive_field_obstruction':{'operator_inequality':'H_field-E_gamma >=0 on the full kinematic field; hence H_physical >= W(H_native tensor I)W*.',
        'held_out_alpha_mass_excess_lower_GeV_exact':str(alpha_gap),
        'any_positive_field_or_positive_kinetic_addition_can_lower_alpha_to_target':False,
        'no_negative_species_dependent_counterterm_inserted':True,
        'polarization_remainder_replaced_by_global_positive_operator_bound':True},
      'complete_records':{'preparation':'Full native uniform fixed-N singlet input and n=0 field service; quartet success, entire matter complement, every nonzero-flux failure and constraint failure retained.',
        'energy_readout':'Spectral measure of the full autonomous operator; all energy complements retained, zero nonselective energy work.',
        'field_preparation_work':'The Gauss-dressed n=0 service has field energy Tr(rho_matter V0)+b, hence field energy above the true vacuum Tr(rho_matter V0)+b-E_gamma in [Tr(rho_matter V0),Tr(rho_matter V0)+b].',
        'zero_charge_vacuum_field_service_only':'For rho supported on all Q_x=0, V0=0 and this contribution lies in [0,b]. This is not the total charged preparation work.',
        'general_full_matter_input':'Retain its actual native charge-square expectation Tr(rho_matter V0); the uniform-quartet values do not replace this full-input expectation.',
        'state_support_not_Hilbert_cutoff':True,'arbitrary_idle_reference_retained':True},
      'complete_infinite_volume_cluster_thresholds':'Not defined on these finite cycles; no asymptotic breakup or halo result is inferred from finite boxes.',
      'physical_nuclear_spectrum_or_common_alpha_attachment_proved':False,
      'new_mass_fit':False,'new_match_with_nature':False,
      'scope':'An exact compact-Gauss extension with the inherited electron conversion, full dynamical cycle field and rigorous native nuclear ground bounds. The positive-extension class cannot cure the alpha mass excess; a Coulomb kernel or a retuned self-energy is not substituted for the full field.'}
