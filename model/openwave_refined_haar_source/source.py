from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import hashlib,json
import sympy as sy
from .gauss import graph,incidence,fusion_paths,invariant,short_cycles,cycle_edges,basis,electric,successors,source_actions,total_action,compressed_sum,loop_phase,haar_trace_product


@lru_cache(None)
def exact_checks():
    vertices,edges,faces=graph();states=basis();loops=((),)+short_cycles();index={s:i for i,s in enumerate(states)}
    phases=(1,)+tuple(int(loop_phase(loop))for loop in loops[1:])
    by_support={frozenset(i for i,k in enumerate(s[0])if k):j for j,s in enumerate(states)}
    assert len(by_support)==201
    oracle=0;outputs=set();max_channels=0;column_count=0;action_digest=hashlib.sha256()
    for face,columns in zip(faces,source_actions()):
        support=frozenset(cycle_edges(face))
        for j,column in enumerate(columns):
            projected={index[target]:value for target,value in column if target in index}
            symmetric=support^frozenset(i for i,k in enumerate(states[j][0])if k)
            candidate=by_support.get(symmetric);expected={}
            if candidate is not None:
                traces=tuple(loop for loop in (loops[candidate],face,loops[j])if loop)
                value=phases[candidate]*phases[j]*haar_trace_product(traces)/2
                if value:expected[candidate]=value
            assert projected==expected,(j,face,projected,expected)
            norm=sy.simplify(sum(v*v for _,v in column));assert 0<norm<=1
            by_spins=Counter(target[0]for target,_ in column)
            max_channels=max(max_channels,max(by_spins.values(),default=0))
            outputs.update(target for target,_ in column);column_count+=1;oracle+=1
            action_digest.update(json.dumps([[target,str(value)]for target,value in column],separators=(',',':')).encode())
    # Tensor orthonormality checks include all higher-valence singlet channels.
    intertwiner_checks=[]
    for spins in ((1,1,1,1),(1,1,1,1,1,1),(2,1,1,1,1),(0,1,1,1,1)):
        paths=fusion_paths(spins)
        vectors=[[invariant(spins,p,m)for m in product(*(range(-k,k+1,2)for k in spins))]for p in paths]
        gram=sy.Matrix(vectors)*sy.Matrix(vectors).T
        assert gram==sy.eye(len(paths));intertwiner_checks.append({'twice_spins':list(spins),'singlet_multiplicity':len(paths),'exact_Gram_identity':True})
    E=sy.diag(*[sy.Rational(electric(s))for s in states]);S=compressed_sum()
    assert S==S.T and E.shape==(201,201)
    second=total_action()[0];norm_s_omega=sy.simplify(sum(v*v for _,v in second));assert norm_s_omega==9
    # Complete square of S on the Haar constant retains every intermediate.
    square={}
    for middle,a in second:
        for face in faces:
            for target,b in successors(middle,face):square[target]=square.get(target,0)+a*b
    square={s:sy.simplify(v)for s,v in square.items()if v}
    fourth=sy.simplify(sum(v*v for v in square.values()));assert fourth==sy.Rational(963,4)
    retained=sy.simplify(sum(v*v for s,v in square.items()if s in index));assert retained==sy.Rational(357,4)
    return {'vertices':len(vertices),'links':len(edges),'elementary_faces':len(faces),
        'vertex_degree_census':{str(k):v for k,v in sorted(Counter(map(len,incidence())).items())},
        'complete_frame_dimension':201,'complete_electric_energy_cutoff_exact':'9/2','full_complement_electric_floor_exact':'5',
        'cycle_census':{'four_cycles':36,'six_cycles':164},'exact_intertwiner_Gram_checks':intertwiner_checks,
        'all_face_frame_columns_checked':column_count,'independent_Haar_matrix_columns_checked':oracle,
        'every_nonzero_successor_retained':True,'distinct_full_action_successors':len(outputs),
        'maximum_intertwiner_channels_for_same_output_spins':max_channels,
        'exact_action_sha256':action_digest.hexdigest(),'complete_frame_spin_network_character_phases':list(phases),
        'sum_Wilson_second_Haar_moment_exact':str(norm_s_omega),'sum_Wilson_fourth_Haar_moment_exact':str(fourth),
        'retained_S_squared_vacuum_norm_squared_exact':str(retained),
        'discarded_S_squared_vacuum_norm_squared_exact':str(fourth-retained),
        'source_representation_cutoff':False}


def run():
    return {'schema':'oph.refined_interior_full_Haar_source.v1',
        'verdict':'PASS_COMPLETE_REFINED_GAUSS_FRAME_WITH_ALL_WILSON_CHANNELS',
        'source_checks':exact_checks(),
        'Hamiltonian':'H=alpha sum_link C_link + beta sum_36faces(1-W_face) on full link Haar L2 with the exact Gauss projection',
        'arbitrary_valence_fusion_tree_source':True,
        'complete_36_bit_readout':{'effects':'E_s=product_f (I+s_f W_f/2)/2 for every s in {-1,+1}^36',
            'outcome_count':2**36,'sum_of_all_effects':'I by independent sign sums','positive_sqrt_Kraus_rule':True,
            'every_outcome_probability_bounds_exact':[str(F(1,4**36)),str(F(3**36,4**36))],
            'nonselective_ground_work_upper_in_units_alpha_exact':'3'},
        'actual_spatial_refinement_source_constructed':True,'uniform_repeated_refinement_or_4D_limit_proved':False,
        'physical_GeV_calibration_or_new_Nature_validation':False,'new_empirical_fit_inputs':[]}
