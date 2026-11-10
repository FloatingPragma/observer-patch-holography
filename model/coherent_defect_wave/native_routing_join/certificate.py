"""Source-pinned common-clock interaction, with complete probability intervals."""
import hashlib
from pathlib import Path
from itertools import product
import sympy as s

from model.particle_common_contract.provenance import import_closure
from ..mobile_flux.routing import object_step
from .model import (State,PROGRAM_ID,encode_lane,native_column,label_word,
                    label_memory,observer,mode,seed,route,native_generator,
                    complete_binary,changed_face,weak_taylor,taylor_error,
                    coin_at_cell,clifford,fixed_native_number)

ROOT=Path(__file__).resolve().parents[3]


def input_hashes():
    here=Path(__file__).parent
    paths={str(p.relative_to(ROOT)) for p in here.iterdir()
           if p.suffix in ('.py','.md') and not p.name.startswith('test_')}
    paths.update(('model/requirements.txt',
                  'model/coherent_defect_wave/mobile_flux/PROOF.md',
                  'model/coherent_defect_wave/mobile_flux/ROUTING_COMPLETION.md',
                  'model/coherent_defect_wave/species/PROOF.md',
                  'model/coherent_defect_wave/spatial/ROUTE.md',
                  'evidence/native_geometric_source_20260925/spec.json',
                  'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
                  'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py'))
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(import_closure(paths)) if (ROOT/name).is_file()}


def experiment():
    count=0
    for label in range(11):
        for occupation in ((),(mode(5,3,0,0),)):
            encoded=encode_lane(label,occupation)
            assert encoded.norm_squared()==1 and observer().gamma(encoded)==encoded
            image=dict(native_column(label,occupation))
            assert s.simplify(sum(s.conjugate(v)*v for v in image.values()))==1
            assert all(int(target!=0)==int(label!=0) and len(n)==len(occupation)
                       and label_memory(target)==label_memory(label) for target,n in image)
            if not occupation:assert image=={(label,()):1}
            count+=1
    inverse_witnesses=((0,(mode(5,3,0,0),)),(1,(mode(7,3,0,0),)),(4,(mode(2,0,0,0),)))
    for label,occupation in inverse_witnesses:
        source=encode_lane(label,occupation)
        assert observer().inverse_step(observer().step(source))==source
        native_column(label,occupation,True)
    # Every collision occupation is retained even with nonempty native registers.
    _,beta=clifford();coin=s.sqrt(3)*s.eye(4)/2-s.I*beta/2
    spectators=((mode(5,3,0,0),),(mode(7,3,0,0),))
    for labels in product(range(11),repeat=2):
        state=State({(labels,spectators):1})
        image=coin_at_cell(state,0,coin)
        assert image.norm_squared()==1
        assert coin_at_cell(image,0,coin.H)==state
        assert all(n==spectators for (_,n),_ in image.amplitudes)
        if all(labels):assert image==state
    spacing=s.pi/12;coupling=3/(2500*s.pi);epsilon=s.Rational(1,10000)
    initial=seed();routed=route(initial,3,spacing)
    assert routed.norm_squared()==1 and fixed_native_number(routed)==1
    assert route(routed,3,spacing,True)==initial
    assert dict(complete_binary(routed,changed_face))[True].norm_squared()==0
    generated=native_generator(routed)
    changed=dict(complete_binary(generated,changed_face))[True]
    assert generated.norm_squared()==s.Rational(283,72)
    assert changed.norm_squared()==s.Rational(1,19200)
    approximate=weak_taylor(routed,epsilon,1)
    delta=taylor_error(initial,epsilon,1)
    assert delta==s.Rational(1,50000000)
    assert approximate.norm_squared()==s.Rational(7200000283,7200000000)
    branches=complete_binary(approximate,changed_face)
    approximate_yes=dict(branches)[True].norm_squared()
    assert approximate_yes==s.Rational(1,1920000000000)
    assert 138**2<19200<139**2
    lower=(epsilon/139-delta)**2;upper=(epsilon/138+delta)**2
    assert epsilon/139>delta and 0<lower<approximate_yes<upper<1
    free=seed(with_native=False);free_routed=route(free,3,spacing)
    assert native_generator(free).norm_squared()==0
    assert native_generator(free_routed).norm_squared()==0
    assert weak_taylor(free_routed,epsilon,3)==free_routed
    cells=tuple(product(range(3),repeat=3));converted={}
    for (labels,natives),value in free_routed.amplitudes:
        occupied=[(lane,label) for lane,label in enumerate(labels) if label]
        assert len(occupied)==1 and not any(natives)
        lane,label=occupied[0]
        converted[cells[lane//2],2*(lane%2)+(label-1)%2,label_word(label)]=value
    expected=object_step({((1,1,1),0,label_word(1)):1},3,spacing)
    assert set(converted)==set(expected)
    assert all(s.simplify(converted[key]-expected[key])==0 for key in expected)
    return {'native_source_encoding_columns_checked':count,
            'literal_native_inverse_witnesses':len(inverse_witnesses),
            'all_columns_checked_by_full_source_reencoding':True,
            'complete_local_routing_columns_with_native_spectators':121,
            'coincident_object_states_retained':100,
            'spacing':'pi/12','coupling':'3/(2500*pi)','weak_duration':'1/10000',
            'cube_side':3,'initial_flux_position':[1,1,1],
            'native_receiver_position':[2,2,2],'native_receiver_lane':0,
            'native_mode':[7,3,0,0],
            'input_norm_squared':'1','native_number':1,'flux_pair_number':1,
            'routed_support':len(routed.amplitudes),'generator_support':len(generated.amplitudes),
            'generator_norm_squared':str(generated.norm_squared()),
            'changed_face_generator_norm_squared':str(changed.norm_squared()),
            'Taylor_order':1,'Taylor_vector_error_upper':str(delta),
            'Taylor_output_norm_squared':str(approximate.norm_squared()),
            'Taylor_output_renormalized':False,
            'Taylor_approximate_projector_weights':[
                {'changed_face':yes,'approximate_squared_norm':str(branch.norm_squared())}
                for yes,branch in branches],
            'exact_unitary_probability_enclosures':[
                {'changed_face':False,'lower':str(1-upper),'upper':str(1-lower)},
                {'changed_face':True,'lower':str(lower),'upper':str(upper)}],
            'complete_binary_instrument':True,'exact_outcome_probabilities_sum_to_one':True,
            'positive_change_lower_bound':str(lower),
            'native_empty_exact_changed_face_probability':'0',
            'coupling_zero_exact_changed_face_probability':'0',
            'free_route_full_output_equal_on_executed_empty_native_input':True,
            'Gauss_preparation':'deterministic code isometry; alternative raw-word heralding retains every sign'}


def certificate():
    return {'schema':'oph.coherent_defect_wave.native_routing_join.v1',
            'program':PROGRAM_ID,
            'shared_tuple':{
                'source':'two routed literal-link lanes with independent 7200-mode native CAR registers per cell',
                'reference':'flat link orbit, empty routing memory and empty native registers at every lane',
                'clock':'J_(a,lambda)=exp(-i*a*lambda*sum_l K_l) R_a; K_l=I-Re(U_F,l) in exact Gauss code',
                'native_U_reads_the_same_link_register_as_routing':True,
                'preparation':'specified normalized code state, coherently prepared by gauge-pair Hadamard and native sign',
                'detector':'complete actual face-syndrome binary PVM with both outcomes',
                'individual_native_number_conserved':True,'total_flux_pair_number_conserved':True,
                'native_registers_held_at_their_spatial_lanes':True},
            'finite_join':experiment(),
            'bounds':{'native_generator':'0<=K<=2N on fixed total native number N',
                      'Taylor_operator_error':'(2*N*abs(a*lambda))**(r+1)/(r+1)!',
                      'normalized_input_required_for_same_vector_error':True},
            'free_restriction':{
                'native_empty_sector_invariant':True,
                'full_finite_route_recovered_on_native_empty_sector':True,
                'one_object_free_frequency':'+-sqrt(|p|^2+4)',
                'rest_frequency_magnitude':2,
                'free_window_error':'t*a*(3*P^2+6*P) for |p_i|<=P',
                'native_empty_self_energy_from_this_coupling':'0',
                'vacuum_particle_production_from_this_clock':'0'},
            'scope':{'finite_common_source_state_clock_join_proved':True,
                     'actual_native_influence_on_routed_flux_proved':True,
                     'complete_source_constraint_and_collision_sectors_retained':True,
                     'weak_native_generator_is_declared_architecture':True,
                     'physical_interacting_relativistic_limit':False,
                     'native_covariant_spatial_dynamics':False,
                     'Lorentz_scalar_interaction_proved':False,
                     'interacting_momentum_window_bound':False,
                     'positive_energy_many_body_vacuum':False,
                     'physical_statistics_or_chiral_matter':False,
                     'scattering_or_SM_identity_or_GeV_mass':False},
            'input_sha256':input_hashes()}
