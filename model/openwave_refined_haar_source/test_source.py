import json
import sympy as sy
from .gauss import graph,basis,fusion_paths,invariant,source_actions,short_cycles,successors
from .source import exact_checks


def test_complete_refined_census_and_full_successor_oracle():
    result=exact_checks()
    assert json.loads(json.dumps(result))==result
    assert result['complete_frame_dimension']==201
    assert result['full_complement_electric_floor_exact']=='5'
    assert result['all_face_frame_columns_checked']==36*201
    assert result['independent_Haar_matrix_columns_checked']==36*201
    assert result['maximum_intertwiner_channels_for_same_output_spins']>1
    assert result['discarded_S_squared_vacuum_norm_squared_exact']=='303/2'


def test_degree_six_multiplicity_is_not_a_single_controller_state():
    paths=fusion_paths((1,)*6)
    assert len(paths)==5
    assert invariant((1,)*6,paths[0],(1,1,1,1,1,1))==0
    assert len(fusion_paths((1,)*5))==0
    assert len(fusion_paths((2,1,1,1,1)))==3


def test_face_products_commute_with_all_intermediate_channels():
    faces=graph()[2];states=basis()
    pairs=((0,1),(0,12),(4,16),(13,25))
    for a,b in pairs:
        state=states[1+(a+b)%36];results=[]
        for first,second in((a,b),(b,a)):
            output={}
            for middle,x in successors(state,faces[first]):
                for target,y in successors(middle,faces[second]):output[target]=output.get(target,0)+x*y
            results.append(output)
        assert all(sy.simplify(results[0].get(k,0)-results[1].get(k,0))==0 for k in results[0].keys()|results[1].keys())
