"""Exact source energies, finite interaction and positive-clock bounds."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import sympy as s

from model.particle_common_contract.provenance import import_closure
from ..background.model import background
from ..gauss import gamma_amplitudes
from .model import (JointObserver,FockObserver,FockState,SparseState,combine,
                    local_action,fock_local_action,dark_wave,dark_pair,
                    dark_plane_wave,energy_moments,lattice,linear_combination,
                    pair_energy_action,mode,MODES)
from ..species.spatial_fock import positive_basis,gamma as pair_gamma
from .dark_collision import certificate as collision_certificate

ROOT=Path(__file__).resolve().parents[3]


def input_hashes():
    local=Path(__file__).parent
    paths={p for p in local.glob('*.py') if not p.name.startswith('test_')}
    paths.update(local/name for name in ('PROOF.md','README.md'))
    paths.update(ROOT/name for name in (
        'model/requirements.txt','model/particle_common_contract/provenance.py',
        'Lean/QFT/PeriodicMomentum.lean','Lean/QFT/SourceHistoryGNSDynamics.lean',
        'evidence/native_geometric_source_20260925/spec.json',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py'))
    names=import_closure({str(p.relative_to(ROOT)) for p in paths})
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(names) if (ROOT/name).is_file()}


def certificate():
    observer=JointObserver();star=observer.law.vertex_stars[5]
    for word in observer.link_words:
        for internal in range(5):
            state=dark_wave(observer,internal,word)
            assert observer.step(state)==combine(state,state,-1,0)
            assert local_action(observer,state)==combine(state,state,2,0)
    state=dark_wave(observer);positive=state.as_dict()
    for key,value in state.amplitudes:
        positive[(key[0]^star,*key[1:])]=value
    positive=SparseState(positive,4)
    assert positive.norm_squared()==1
    assert gamma_amplitudes(positive.as_dict())==positive.as_dict()
    assert local_action(observer,positive)==combine(positive,positive,2,0)
    zero=background(observer,0)
    assert local_action(observer,zero).norm_squared()==0
    fock=FockObserver();pair=dark_pair(fock);values=pair.as_dict()
    for (word,wedge),value in pair.amplitudes:
        values[word^star,wedge]=value
    pair=FockState(values,8)
    assert pair.norm_squared()==1 and fock.gamma(pair)==pair
    assert fock.step(pair)==pair and fock_local_action(fock,pair).norm_squared()==0
    two=positive_basis((1,3),(mode(5,3,0,0),MODES+mode(6,0,0,0)))
    two_image=pair_energy_action(two);lookup=two.as_dict()
    two_mean=s.expand(sum(s.conjugate(lookup.get(key,0))*value
                         for key,value in two_image.amplitudes)/two.norm_denominator)
    two_second=two_image.norm_squared()
    assert two_mean==s.Rational(89,30)+s.pi/2
    assert s.expand(two_second-(s.Rational(3371,360)+89*s.pi/30+3*s.pi**2/8))==0
    assert all(pair_gamma(two_image,c)==two_image for c in (0,1))
    seed=lattice.packet();moments=energy_moments(seed,4)
    expected={'mean':s.Rational(59,30)+3*s.pi/4,
              'second':s.Rational(283,72)+59*s.pi/20+3*s.pi**2/4,
              'variance':s.Rational(113,1800)+3*s.pi**2/16}
    assert all(s.expand(moments[key]-value)==0 for key,value in expected.items())
    image=moments['A_image']
    moved=lattice.State({key:value for key,value in image.amplitudes if key[2]},image.norm_denominator)
    q=moved.norm_squared();assert q==F(1,1800)
    t=F(1,10000);bound=F(80,7);remainder=t*t*bound*bound/2
    lower=(t/43-remainder)**2;upper=(t/42+remainder)**2
    assert t/43>remainder and 0<lower<t*t*q<upper<1
    planes=[]
    for momentum in ((0,0,0),(1,0,0),(1,1,0),(1,1,1),(2,0,0),(2,2,2)):
        plane=dark_plane_wave(4,momentum);result=energy_moments(plane,4)
        energy=2+s.pi/2*sum(1-s.cos(s.pi*k/2) for k in momentum)
        assert result['variance']==0 and s.expand(result['mean']-energy)==0
        a,b=result['A_image'],result['hopping_coefficient_image']
        assert a==linear_combination((plane,),(2,))
        lattice_energy=sum(1-s.cos(s.pi*k/2) for k in momentum)/2
        assert b==linear_combination((plane,),(F(lattice_energy),))
        planes.append({'momentum_in_pi_over_two_units':list(momentum),'energy':str(energy),
                       'energy_variance':'0','occupied_wave_support':len(plane.amplitudes)})
    return {'schema':'oph.coherent_defect_wave.positive_energy.v1',
            'program':'oph.coherent_defect_wave.positive_continuous_clock.v1',
            'input_sha256':input_hashes(),
            'source_energy_law':'H=sum_x[I-(U_F,x+U_F,x*)/2]+(pi/4)sum_edges dGamma(I-S_xy)',
            'continuous_evolution':'exp(-itH)',
            'clock_equals_discrete_seven_layer_tick':False,
            'Gauss_constraints_and_total_number_preserved':True,
            'positive_Gauss_empty_wave_ground_energy':'0',
            'regional_topological_charge_per_cell':'Q_x=I on the admitted link sector',
            'intercell_transport_of_regional_topological_charge':False,
            'wave_energy_bands_identified_as_traveling_protected_fluxes':False,
            'ground_selected_unique':False,
            'fixed_number_operator_bound':'0 <= H <= (2+3*pi)*N',
            'exact_one_wave_operator_norm_on_even_periodic_cubes':'2+3*pi',
            'bound_independent_of_number_of_cells_at_fixed_N':True,
            'dark_U_minus_one_witnesses_checked':50,
            'dark_onsite_energy':'2','uniform_triplet_onsite_energy':'0',
            'positive_Gauss_dark_state_checked':True,
            'same_cell_two_dark_wave_onsite_energy':'0',
            'separated_two_dark_wave_onsite_energy':'4',
            'onsite_energy_is_second_quantized_one_body_energy':False,
            'two_cell_two_wave_energy':{'mean':str(two_mean),'second_moment':str(two_second),
                                       'H_image_support':len(two_image.amplitudes),
                                       'both_Gauss_constraints_preserved':True,
                                       'input_positive_Gauss_preparation_success':'1/4',
                                       'moments_conditional_on_that_preparation':True},
            'cubic_source_packet':{'length':4,'support_cells':8,
                                   'mode':[7,3,0,0],'source': 'cubic_readout.lattice.packet',
                                   'moments':{k:str(moments[k]) for k in expected},
                                   'A_image_support':len(image.amplitudes),
                                   'hopping_coefficient_support':len(moments['hopping_coefficient_image'].amplitudes)},
            'finite_time_flux_measurement':{
                'initial_changed_flux_probability':'0','generator_changed_flux_norm_squared':str(q),
                'time':str(t),'operator_norm_upper_bound':str(bound),
                'unitary_first_order_remainder_upper_bound':str(remainder),
                'yes_probability_lower_bound':str(lower),'yes_probability_upper_bound':str(upper),
                'no_probability_interval':[str(1-upper),str(1-lower)],
                'all_binary_outcomes_retained':True,
                'probability_enclosure_from_executed_exact_generator':True,
                'full_time_histogram_sampled':False,
                'conditioned_on_normalized_positive_Gauss_preparation':True},
            'exact_dark_band':'2+(pi/2)*sum_i(1-cos(k_i))',
            'exact_dark_plane_waves':planes,
            'finite_dark_collision':collision_certificate(),
            'energy_estimator':{'onsite_Kraus':'(I+U_F)/2, (I-U_F)/2',
                                'onsite_record_values':['0','2'],
                                'edge_record_values':'(pi/2)*antisymmetric_mode_occupation',
                                'noncommuting_term_measurements_form_joint_energy_PVM':False},
            'scope':{'local_positive_energy_convention_specified':True,
                     'finite_occupation_interaction':True,
                     'physical_clock_in_seconds':False,'relativistic_dispersion_proved':False,
                     'physical_rest_mass_identified':False,'GeV_prediction':False}}
