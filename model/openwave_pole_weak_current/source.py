"""Actual native-current matrix elements between isolated source poles."""
from fractions import Fraction as F
from functools import lru_cache
import sympy as sy
from model.openwave_baryon_poles import source as poles
from model.openwave_baryon_dispersion.exact_native import complete
from model.openwave_hadron_refinement_flow.analytic import constants
from model.openwave_hadron_refinement_flow.parameters import degree_coefficients
from model.openwave_hadron_refinement_flow.certified_rest import full,sqrt_upper
from model.openwave_capacity_stress.source import sqrt_interval,enclosure
from .native import contraction


def transfer_error(neutron,proton,norm):
    pn,pp=neutron['pole'],proton['pole']
    qn,qp=(F(x['full_Q_probability_upper'])for x in(pn,pp))
    an,ap=(F(x['P_ground_sine_upper'])+F(x['prepared_trial_sine_upper'])for x in(pn,pp))
    return F(norm)*(sqrt_upper(2)*(an+ap)+qn+qp+sqrt_upper(qn*qp))


def absolute_interval(numerator,normalization,error,norm):
    square=poles.real_fraction(sy.expand(sy.conjugate(numerator)*numerator))/normalization
    lo,hi=sqrt_interval(square,24)
    return square,(max(F(0),lo-error),min(F(norm),hi+error))


@lru_cache(None)
def work_bound():
    A,G,_=complete();g=list(G.diagonal())
    floor=min(F(A[i,i])/F(g[i])-sum(abs(F(A[i,j]))*sqrt_upper(F(1,F(g[i])*F(g[j])))for j in range(100)if i!=j)for i in range(100))
    assert floor>0
    _,_,zeta=poles.bands.repair_matrix(False)
    _,_,_,p,_,unit=constants(*degree_coefficients())
    spin=F(27,4)*(abs(p['Js']+3*p['spin_gradient'])+abs(p['Jw']-p['eta'])+2*abs(p['eta'])+2*abs(p['spin_gradient']))
    R=F(3,2)*abs(p['delta'])+spin+4*zeta/floor
    return {'noncommuting_bounded_source_remainder_norm_upper':str(R),
        'complete_instrument_mean_work_absolute_upper_source':str(2*R),
        'complete_instrument_mean_work_absolute_upper_GeV':str(2*R*unit),
        'unbounded_colour_and_number_potential_commutes_with_instrument':True,
        'argument':'The complete light three-quark source equals a commuting scalar/colour/hopping part plus a bounded spin-isospin, charge-split and correlated-kinetic remainder. The dual complete instrument is unital and contracts its norm.'}


def point(cosine,sine):
    n=poles.phase(False,cosine,sine);p=poles.phase(True,cosine,sine)
    vn=poles.vector(n['trial_source_coordinates']);vp=poles.vector(p['trial_source_coordinates'])
    nn=F(n['trial_Gram_norm_squared']);np=F(p['trial_Gram_norm_squared']);data=contraction();rows={}
    for kind,norm in(('vector',1),('axial',3)):
        numerator=sy.expand((vp.conjugate().T*data[kind+'_numerator']*vn)[0])
        error=transfer_error(n,p,norm)
        squared,interval=absolute_interval(numerator,nn*np,error,norm)
        origin='Complete native-current contraction with exact compact/trial/tail transfer error.'
        if cosine==1 and sine==0:
            if kind=='vector':interval=(F(1),F(1));origin='Exact isospin raising intertwines the full rest neutron and proton sources.'
            else:interval=tuple(full(*degree_coefficients())['full_source_axial_interval']);origin='Exact same-state rest transfer with the parent native axial-variance bound.'
        assert interval[0]>0
        rows[kind]={'actual_source_pole_amplitude_magnitude_interval':list(map(str,interval)),
            'actual_source_pole_amplitude_display':enclosure(interval),
            'normalized_rational_trial_amplitude_squared_exact':str(squared),
            'general_compact_to_full_amplitude_error_upper':str(error),
            'current_success_then_proton_pole_probability_interval':[str(interval[0]**2/9),str(interval[1]**2/9)],
            'original_trial_input_neutron_pole_then_current_then_proton_pole_probability_lower':str(F(n['pole']['prepared_trial_pole_weight_lower'])*interval[0]**2/45),
            'certificate':origin}
    return {'cosine':str(cosine),'sine':str(sine),'neutron_source_pole_GeV_interval':n['source_pole_GeV_interval'],
        'proton_source_pole_GeV_interval':p['source_pole_GeV_interval'],'currents':rows}


def run():
    data=contraction();rows=[point(c,s)for c,s in poles.PHASES]
    _,_,_,p,_,unit=constants(*degree_coefficients())
    return {'schema':'oph.native_pole_weak_current.v1','verdict':'PASS_NATIVE_WEAK_CURRENTS_ON_PROPAGATING_STRONG_SOURCE_POLES',
        'momentum_records':rows,
        'exact_native_current':{'all_compact_triquark_shapes':data['complete_native_compact_shapes'],
            'complete_fixed_spin_charge_shapes':data['neutron_fixed_weight_shapes'],
            'vector_five_mode_Gram_numerator':[[str(x)for x in row]for row in data['vector_numerator'].tolist()],
            'axial_five_mode_Gram_numerator':[[str(x)for x in row]for row in data['axial_numerator'].tolist()],
            'complete_fixed_weight_axial_square_matrix':[[str(x)for x in row]for row in data['axial_full_fixed_weight_square'].tolist()],
            'complete_fixed_weight_axial_square_spectrum_contains_only':['1','4'],
            'current_flavour_or_spin_output_discarded':False,'floating_rounding_supplies_current_entries':False},
        'complete_native_instrument':{'success':'K_success=J/3 on the complete light three-quark source, with a charge-minus-one reference shift; zero success on its full Fock complement.',
            'failure':'K_failure=sqrt(I-J*J/9), with identity failure on the full other-sector complement.',
            'reference':'Two-sided integer charge reference with declared zero Hamiltonian; every final reference charge is retained.',
            'source_charge_change_on_success_exact':'1','reference_charge_change_on_success_exact':'-1','total_charge_change_exact':'0',
            'terminal_outcomes':'proton source pole, every excited/current-produced source sector and preparation/current failure',
            'terminal_spectral_readout_work_exact':'0','arbitrary_idle_external_quantum_reference':True,
            'rest_vector_success_probability_exact':'1/9','rest_vector_failure_probability_exact':'8/9',
            'rest_vector_mean_intervention_work_GeV_exact':str(-p['delta']*unit/9),
            'full_source_work_bound':work_bound()},
        'new_empirical_current_or_mass_input':False,
        'ancestry':'Inherited strong mass/axial calibrations, electron GeV anchor and declared positive propagation rule. The native current is unchanged.',
        'autonomous_weak_decay_or_neutron_lifetime_proved':False,
        'physical_3plus1_current_or_particle_dictionary_proved':False,
        'scope':'Actual native weak-current amplitudes between isolated fixed-spacing source poles. The current intervention retains its work and reference; no autonomous beta decay or physical lifetime is inferred.'}
