"""Exact moving-native records, source spin extension and scaling boundaries."""
import hashlib
from pathlib import Path
import sympy as s

from model.particle_common_contract.provenance import import_closure
from ..native_routing_join.model import native_column as primitive_column
from ..mobile_flux.model import gamma as literal_gamma
from . import model as m
from . import spinful as f

ROOT=Path(__file__).resolve().parents[3]


def input_hashes():
    here=Path(__file__).parent
    paths={str(p.relative_to(ROOT)) for p in here.iterdir()
           if p.suffix in ('.py','.md') and not p.name.startswith('test_')}
    paths.update(('model/requirements.txt',
        'model/coherent_defect_wave/native_routing_join/PROOF.md',
        'model/coherent_defect_wave/mobile_flux/ROUTING_COMPLETION.md',
        'model/coherent_defect_wave/species/PROOF.md',
        'model/coherent_defect_wave/spatial/ROUTE.md',
        'evidence/native_geometric_source_20260925/spec.json',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py'))
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(import_closure(paths)) if (ROOT/name).is_file()}


def matching_experiment():
    i=m.mode(5,3,0,0);j=m.mode(6,0,0,0)
    source=m.State({((1,0),((i,),(j,))):1})
    output=m.edge(source,0,1)
    assert source.norm_squared()==output.norm_squared()==1
    assert m.edge(output,0,1,True)==source
    encoded=m.encode_pair(source);literal=m.literal_pair_edge(encoded)
    assert literal==m.encode_pair(output)
    assert all(literal_gamma(literal,lane)==literal for lane in (0,1))
    occupation=[m.State({key:v for key,v in output.amplitudes if len(key[1][0])==n}).norm_squared()
                for n in range(3)]
    assert occupation==[s.Rational(1,4),s.Rational(1,2),s.Rational(1,4)]
    identical=m.State({((1,0),((j,),(j,))):1})
    assert m.edge(identical,0,1)==m.State({((1,0),((j,),(j,))):-s.I})
    state,receiver,sender=m.moving_seed();routed=m.route(state,4,s.pi/12)
    moved=m.native_matching(routed,4,0)
    assert moved.norm_squared()==1 and m.native_matching(moved,4,0,True)==routed
    generator=m.native_generator(moved)
    assert generator.norm_squared()==s.Rational(283,72)
    records=dict(m.complete_records(generator,receiver))
    assert records[True,True].norm_squared()==s.Rational(1,38400)
    assert records[False,True].norm_squared()==0
    epsilon=s.Rational(1,10000);delta=s.Rational(1,50000000)
    assert 195**2<38400<196**2 and epsilon/196>delta
    lower=(epsilon/196-delta)**2;upper=(epsilon/195+delta)**2
    approximate=m.weak_taylor(moved,epsilon,1)
    assert approximate.norm_squared()==s.Rational(7200000283,7200000000)
    return {'program':m.PROGRAM_ID,
            'clock':'six-state scheduler: exp(-ia lambda K) T_phase R_a; phase advances modulo6',
            'native_one_body_modes_per_lane':7200,
            'literal_Gauss_intertwiner_checked':True,
            'full_CAR_pair_number_records':[str(x) for x in occupation],
            'identical_mode_pair_phase':'-I',
            'native_sender_position':[3,2,2],'native_receiver_position':[2,2,2],
            'source_cube_side':4,'phase':0,'next_phase':1,
            'native_arrival_probability':'1/2',
            'routed_and_moving_support':[len(routed.amplitudes),len(moved.amplitudes)],
            'generator_support':len(generator.amplitudes),
            'changed_face_generator_norm_squared':'1/38400',
            'weak_duration':str(epsilon),'Taylor_vector_error':str(delta),
            'Taylor_output_norm_squared':str(approximate.norm_squared()),
            'Taylor_output_renormalized':False,
            'exact_joint_records':[
                {'arrived':False,'changed':False,'lower':'1/2','upper':'1/2'},
                {'arrived':False,'changed':True,'lower':'0','upper':'0'},
                {'arrived':True,'changed':False,'lower':str(s.Rational(1,2)-upper),'upper':str(s.Rational(1,2)-lower)},
                {'arrived':True,'changed':True,'lower':str(lower),'upper':str(upper)}],
            'primitive_matching_relativistic_native_limit':False}


def spinful_experiment():
    primitive=f.mode(7,3,0,0)
    for spin in range(4):
        got=dict(f.native_column(1,(4*primitive+spin,)))
        expected={(label,(4*n[0]+spin,)):v for (label,n),v in primitive_column(1,(primitive,))}
        assert got==expected
        assert f.gamma(f.encode(1,(4*primitive+spin,)))==f.encode(1,(4*primitive+spin,))
    pair=f.SpinFock({(1,(4*primitive,4*primitive+1)):1})
    assert f.controls(pair)==pair
    initial,target,receiver=f.focused_preparation()
    assert f.norm(initial)==1 and len(initial)==8
    assert f.complete_records(initial,receiver)[False,False]==1
    assert f.native_route(initial,4)==target
    routed=f.flux_route(initial,4,s.pi/12)
    moved=f.native_route(routed,4)
    assert f.norm(moved)==1 and f.complete_records(moved,receiver)[True,False]==1
    assert f.norm({key:value for key,value in routed.items()
                   if any(label and any(index//f.EXT_MODES==lane for index in key[1])
                          for lane,label in enumerate(key[0]))})==0
    # Full exterior two-native propagation includes Pauli signs and is reversible.
    labels=next(iter(target))[0]
    pair_state={(labels,(4*f.mode(2,0,0,0),4*f.mode(2,0,0,0)+1)):s.Integer(1)}
    pair_moved=f.native_route(pair_state,4)
    assert f.norm(pair_moved)==1 and f.native_route(pair_moved,4,True)==pair_state
    k=f.generator(moved);v=f.contact(moved)
    assert f.norm(k)==s.Rational(283,72) and f.norm(v)==s.Rational(1,9600)
    assert f.complete_records(k,receiver)[True,True]==s.Rational(1,19200)
    assert f.complete_records(v,receiver)[True,True]==s.Rational(1,19200)
    epsilon=s.Rational(1,10000);delta_k=2*epsilon**2;delta_v=8*epsilon**2
    assert 138**2<19200<139**2 and epsilon/139>delta_v
    lower=(epsilon/139-delta_v)**2;upper=(epsilon/138+delta_v)**2
    approx=f.add((moved,v),(1,-s.I*epsilon))
    assert f.norm(approx)==1+epsilon**2/s.Integer(9600)
    return {'program':f.PROGRAM_ID,
            'added_spin_memory_dimension':4,'native_one_body_modes_per_lane':f.EXT_MODES,
            'native_source_law':'exterior(W_g tensor I4), then endpoint parity summed over allspin components on one sharedg',
            'complete_same_link_spin_columns':4,
            'two_spin_occupied_parity_control_checked':True,
            'native_kinetic_law':'Gamma(three Clifford shifts), mass angle zero',
            'flux_kinetic_law':'complete paired-flux route with mass angle2a',
            'clock_K':'exp(-ia lambda K) Gamma(D_a^0) R_a',
            'clock_V':'exp(-ia lambda (K-K0)) Gamma(D_a^0) R_a',
            'K0_definition':'actual flat-orbit native generator extended identically over flux labels',
            'flat_background_subtraction_is_separate_declared_law':True,
            'subtraction_is_scalar_energy_shift_or_Wick_ordering':False,
            'incoming_native_packet_support':len(initial),
            'initial_native_receiver_probability':'0','propagated_native_receiver_probability':'1',
            'native_propagation_off_changed_face_probability':'0',
            'exact_two_native_exterior_inverse_checked':True,
            'complete_two_native_route_support':len(pair_moved),
            'joint_kinetic_output_support':len(moved),
            'K_image_support':len(k),'V_image_support':len(v),
            'K_norm_squared':'283/72','V_norm_squared':'1/9600',
            'changed_face_generator_norm_squared':'1/19200',
            'weak_duration':str(epsilon),'K_Taylor_vector_error':str(delta_k),
            'V_Taylor_vector_error':str(delta_v),
            'V_Taylor_output_norm_squared':str(f.norm(approx)),
            'Taylor_output_renormalized':False,
            'exact_V_clock_joint_records':[
                {'arrived':False,'changed':False,'lower':'0','upper':'0'},
                {'arrived':False,'changed':True,'lower':'0','upper':'0'},
                {'arrived':True,'changed':False,'lower':str(1-upper),'upper':str(1-lower)},
                {'arrived':True,'changed':True,'lower':str(lower),'upper':str(upper)}],
            'positive_changed_face_lower':str(lower),
            'free_native_frequency':'+-abs(p_native)',
            'free_flux_frequency':'+-sqrt(abs(p_flux)^2+4)',
            'combined_free_window_error':'t*a*(3*P_flux^2+6*P_flux+3*P_native^2)',
            'unchanged_K_bulk_frequency':'lambda*kappa_j +-abs(p_native)',
            'primitive_internal_labels_identified_as_added_spin':False}


def exact_scaling_checks():
    q=s.Symbol('q',real=True)
    X=s.Matrix([[0,1],[1,0]]);Xq=s.Matrix([[0,s.exp(-s.I*q)],[s.exp(s.I*q),0]])
    T=(m.A*s.eye(2)+m.B*Xq)*(m.A*s.eye(2)+m.B*X)
    assert s.simplify(s.expand_complex(s.trace(T))-s.I*(s.cos(q)-1))==0
    assert s.simplify(T.det())==-1
    alpha,beta=f.clifford();P=s.diag(1,1,0,0);residuals=[]
    for A in alpha:
        G=s.kronecker_product((A*P+P*A)/2,s.eye(4))+s.kronecker_product(P,A)
        assert s.trace(G.H*G)==16
        B=s.kronecker_product(A,s.eye(4))+s.kronecker_product(s.eye(4),A)
        scalar=s.kronecker_product(beta,beta)
        assert B*scalar+scalar*B==s.zeros(16)
        residuals.append(str(s.trace(G.H*G)))
    return {'primitive_matching_trace':'-2*i*sin(q/2)^2','primitive_matching_determinant':'-1',
            'primitive_low_phase':'arcsin(sin(q/2)^2)',
            'common_clock_native_kinetic_bound':'frequency<=a*abs(p)^2/3 for |a*p_i|<=1/2',
            'primitive_native_kinetic_freezes_in_a_clock':True,
            'contact_domain':'one flux and one native, fixed physical torus, preserved nativeM-mode window',
            'contact_operator_bound':'V=CVC, norm(V)<=4',
            'exact_discrete_Duhamel_bound':'4*a*abs(lambda)*sum_j norm(C*psi_free,j)',
            'Dirac_contact_state_difference_bound':'4*abs(lambda(a))*t*sqrt(M/ell^3)*a^(3/2)',
            'two_cell_matching_extra_factor':'sqrt(8)',
            'interacting_state_replaced_by_free_probability':False,
            'bounded_lambda_smooth_fixed_particle_interaction_survives':False,
            'dense_sea_or_strong_scaling_excluded_by_this_bound':False,
            'standard_spinor_density_vertex_boost_residual_squared':residuals,
            'scalar_beta_tensor_beta_residual':'0',
            'density_vertex_proved_Lorentz_scalar':False}


def certificate():
    return {'schema':'oph.coherent_defect_wave.propagating_native_join.v1',
            'primitive_matching':matching_experiment(),
            'spinful_common_clock':spinful_experiment(),
            'scaling_and_vertex':exact_scaling_checks(),
            'scope':{'complete_native_spatial_transport':True,
                     'controlled_free_native_and_flux_Dirac_windows':True,
                     'complete_finite_shared_link_interaction':True,
                     'interacting_relativistic_limit':False,
                     'positive_energy_relativistic_many_body_reference':False,
                     'source_mass_renormalization_or_SM_identity_or_GeV_prediction':False},
            'input_sha256':input_hashes()}
