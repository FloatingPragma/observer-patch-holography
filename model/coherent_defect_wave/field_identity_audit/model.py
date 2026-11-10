"""Exact source-module attachment and controls against assigning face labels charges."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations,product
import sympy as s

from ..model import SparseState
from ..mobile_flux.model import observer,STAR
from ..mobile_flux.routing import odd_words
from ..native_routing_join.model import native_column
from ..species.fock import FockState,mode,coordinates,wedge_insert
from ..physical_attachment.chiral import representation


def dark_orbital(word,internal):
    return SparseState({(word,2,internal,0,0):F(1),
                        (word,2,internal,0,1):F(-1)},2)


def exterior_word(word,subsets,images=None):
    """Actual exterior polynomial of coherent one-body source columns."""
    terms={():s.Integer(1)}
    for internal in subsets:
        orbital=(images[internal] if images is not None else dark_orbital(word,internal))
        following={}
        for wedge,coefficient in terms.items():
            for (g,p,j,r,a),value in orbital.amplitudes:
                assert g==word
                target,sign=wedge_insert(wedge,mode(p,j,r,a))
                if sign:following[target]=following.get(target,0)+sign*coefficient*value
        terms={key:s.expand(value) for key,value in following.items() if value!=0}
    return FockState({(word,wedge):value for wedge,value in terms.items()},2**len(subsets))


def dark_code(word,subset):
    first=exterior_word(word,subset)
    out=first.as_dict()
    for (g,wedge),value in first.amplitudes:out[g^STAR,wedge]=value
    return FockState(out,2*first.norm_denominator)


def current_action(state,current):
    """dGamma(I_port,repair,address tensor current) on the literal native CAR modes."""
    output={}
    for (word,wedge),coefficient in state.amplitudes:
        for ordinal,index in enumerate(wedge):
            p,j,r,a=coordinates(index)
            rest=wedge[:ordinal]+wedge[ordinal+1:]
            for target_internal in range(5):
                value=current[target_internal,j]
                target=mode(p,target_internal,r,a)
                if value==0 or target in rest:continue
                sign=(-1)**(ordinal+sum(i<target for i in rest))
                key=word,tuple(sorted((*rest,target)))
                output[key]=s.expand(output.get(key,0)+sign*value*coefficient)
    return FockState({key:s.expand(value) for key,value in output.items() if value!=0},state.norm_denominator)


@lru_cache(maxsize=1)
def dark_module():
    source=observer();subsets,g5,currents,_=representation()
    assert source.base.addresses[:2]==((0,1),(1,0))
    assert all(2 not in source.base.law.edges[e] for e in source.base.law.interior)
    body_checks=0;occupation_checks=0
    for word in source.base.link_words:
        images={j:source.base.wave_step(dark_orbital(word,j)) for j in range(5)}
        for j,image in images.items():
            original=dark_orbital(word,j)
            assert image==SparseState({key:-value for key,value in original.amplitudes},2)
            body_checks+=1
        for subset in subsets:
            original=exterior_word(word,subset)
            actual=source.controls(exterior_word(word,subset,images))
            assert actual==original and original.norm_squared()==1
            assert source.controls(original,True)==original
            occupation_checks+=1
    code=[dark_code(1,subset) for subset in subsets]
    assert all(source.gamma(state)==state and state.norm_squared()==1 for state in code)
    assert all(not set(code[i].as_dict())&set(code[j].as_dict())
               for i in range(15) for j in range(i))
    current_checks=0
    for primitive,represented in zip(g5,currents,strict=True):
        for column,state in enumerate(code):
            expected={}
            for row,coefficient in enumerate(represented[:,column]):
                if coefficient==0:continue
                assert code[row].norm_denominator==state.norm_denominator
                for key,value in code[row].amplitudes:
                    expected[key]=s.expand(expected.get(key,0)+coefficient*value)
            actual=current_action(state,primitive)
            assert actual==FockState(expected,state.norm_denominator)
            current_checks+=1
    vacuum=FockState({(1,()):1,(1^STAR,()):1},2)
    assert all(current_action(vacuum,t).norm_squared()==0 for t in g5)
    assert source.step(vacuum)==source.inverse_step(vacuum)==vacuum
    for subset,state in zip(subsets,code,strict=True):
        expected=FockState({key:len(subset)*value for key,value in state.amplitudes},state.norm_denominator)
        assert current_action(state,s.eye(5))==expected
    degree_sums=[sum(len(subset)**power for subset in subsets) for power in (1,3)]
    assert degree_sums==[40,400]
    y=currents[-1]
    assert y.det()!=0
    return {'one_body_source_columns_checked':body_checks,
            'factorized_full_exterior_and_actual_control_columns_checked':occupation_checks,
            'actual_native_current_intertwiner_columns_checked':current_checks,
            'native_mode':'port2, repair0, normalized address0 minus address1',
            'native_occupation_degrees':[2,4],'internal_dimension':15,
            'twelve_current_algebra':'su(3)+su(2)+u(1)_Y inside native u(3)+u(2)',
            'extra_native_degree_current_columns_checked':15,
            'extra_degree_if_gauged_as_all_left_Weyl_anomaly_sums':{'linear':degree_sums[0],'cubic':degree_sums[1]},
            'full_even_dark_module_dimension_including_vacuum':16,
            'omitted_even_vacuum_current_representation':'one gauge singlet',
            'literal_even_native_occupation_module':True,
            'native_U_on_module':'+I','native_K_on_module':'0',
            'module_is_Gauss_positive':True,
            'hypercharges':[str(q) for q in y.diagonal()],
            'basis_subsets':[list(subset) for subset in subsets],
            'hypercharge_operator_invertible':True,
            'empty_native_flux_code_actual_current_representation':'trivial',
            'current_intertwiner_from_empty_native_flux_to_R15_dimension':0,
            'module_prepared_or_bound_by_native_idle_law':False}


@lru_cache(maxsize=1)
def face_charge_obstruction():
    """Actual K transitions obstruct assigning a 3+2 charge split to five face positions."""
    native=(mode(5,3,0,0),);edges={};images={}
    for face in range(5):
        label=1+2*face;out=defaultdict(F);out[label,native]=F(1)
        for inverse in (False,True):
            for key,value in native_column(label,native,inverse):out[key]-=F(value)/2
        images[face]={key:value for key,value in out.items() if value}
        for (target,occupation),value in images[face].items():
            other=(target-1)//2
            if other!=face:
                edges.setdefault(tuple(sorted((face,other))),
                                 {'input_face':face,'output_face':other,
                                  'output_native_coordinates':[list(coordinates(i)) for i in occupation],
                                  'coefficient':str(value)})
    incidence=s.zeros(len(edges),5)
    for row,(a,b) in enumerate(sorted(edges)):incidence[row,a]=1;incidence[row,b]=-1
    assert incidence.rank()==4
    native_blocks={}
    for face,image in images.items():
        for (label,occupation),value in image.items():
            native_blocks.setdefault(occupation,s.zeros(5))[(label-1)//2,face]+=value
    unique={tuple(block):block for block in native_blocks.values()}
    columns=[]
    for i,j in product(range(5),repeat=2):
        unit=s.zeros(5);unit[i,j]=1
        columns.append(s.Matrix.vstack(*[(block*unit-unit*block).reshape(25,1)
                                          for block in unique.values()]))
    pure_face_rank=s.Matrix.hstack(*columns).rank()
    assert pure_face_rank==24
    controls=[]
    for triplet in combinations(range(5),3):
        q=[2 if i in triplet else -3 for i in range(5)]
        norm=sum((value*(q[(label-1)//2]-q[0]))**2
                 for (label,_),value in images[0].items())
        assert norm>0
        controls.append({'triplet_face_indices':list(triplet),'actual_K_charge_commutator_norm_squared':str(norm)})
    return {'source_native_input':[5,3,0,0],'face_transition_graph_edges':[list(edge) for edge in sorted(edges)],
            'nonzero_full_source_matrix_element_witnesses':[edges[edge] for edge in sorted(edges)],
            'diagonal_face_current_commutant_dimension':1,
            'actual_native_output_blocks':len(native_blocks),
            'distinct_face_matrices_in_current_test':len(unique),
            'full_pure_face_current_constraint_rank':pure_face_rank,
            'full_pure_face_current_commutant_dimension':1,
            'all_three_plus_two_face_splits_fail':controls,
            'actual_native_C5_internal_current_identified_with_face_positions':False}


def one_site_current(vector,current):
    output={}
    for (left,right),value in vector.items():
        for cell,label in enumerate((left,right)):
            if label==0:continue
            for target in range(1,16):
                coefficient=current[target-1,label-1]
                if coefficient:
                    pair=(target,right) if cell==0 else (left,target)
                    output[pair]=s.expand(output.get(pair,0)+coefficient*value)
    return {key:value for key,value in output.items() if value!=0}


def graded_swap(vector):
    return {(right,left):(-value if left and right else value)
            for (left,right),value in vector.items()}


def bound_transport():
    """Complete two-slot graded exchange on vacuum plus the actual native R15 cloud."""
    _,_,currents,_=representation();checks=0
    for pair in product(range(16),repeat=2):
        state={pair:s.Integer(1)}
        assert graded_swap(graded_swap(state))==state
        for current in currents:
            assert graded_swap(one_site_current(state,current))==one_site_current(graded_swap(state),current)
            checks+=1
    a,b=(1-s.I)/2,(1+s.I)/2
    assert s.simplify(a*s.conjugate(a)+b*s.conjugate(b))==1
    assert s.simplify(a*s.conjugate(b)+b*s.conjugate(a))==0
    initial={(1,0):s.Integer(1)}
    output={key:a*value for key,value in initial.items()}
    for key,value in graded_swap(initial).items():output[key]=output.get(key,0)+b*value
    probabilities=[s.simplify(s.conjugate(output.get(pair,0))*output.get(pair,0))
                   for pair in product(range(16),repeat=2)]
    assert probabilities[16]==probabilities[1]==s.Rational(1,2)
    assert sum(probabilities)==1
    # A single orbital slot has an occupancy constraint. These projected
    # creation matrices are not falsely called the full fifteen-mode CAR.
    creation=s.zeros(16);creation[1,0]=1
    assert creation*creation.H+creation.H*creation!=s.eye(16)
    assert s.simplify(a-b)==-s.I
    return {'per_slot_internal_code_dimension':16,'complete_two_slot_dimension':256,
            'total_native_current_commutator_columns_checked':checks,
            'swap_is_complete_unitary_and_preserves_total_native_current':True,
            'native_cloud_moves_with_flux_object':True,
            'occupied_flux_grade':'odd declared flux grade times actual even native cloud',
            'two_occupied_slot_exchange_sign':-1,'two_equal_occupied_slots_partial_exchange_phase':'-I',
            'complete_256_outcome_probabilities':list(map(str,probabilities)),
            'physical_preparation':'one bound cloud in the left slot; right slot flat and native-empty',
            'full_unconstrained_composite_CAR_derived_from_two_slots':False,
            'native_idle_law_generates_the_binding_or_transport':False}


def antisymmetric_clouds():
    """Actual sum-of-cloud currents on an antisymmetrized two-object encoding."""
    _,_,currents,_=representation();pairs=tuple(combinations(range(15),2));checks=0
    def encode(pair):
        a,b=pair
        return {(a+1,b+1):1/s.sqrt(2),(b+1,a+1):-1/s.sqrt(2)}
    for pair in pairs:
        source=encode(pair)
        assert sum(s.conjugate(v)*v for v in source.values())==1
        for current in currents:
            target={}
            for ordinal,j in enumerate(pair):
                rest=pair[:ordinal]+pair[ordinal+1:]
                for i in range(15):
                    coefficient=current[i,j]
                    if coefficient==0 or i in rest:continue
                    sign=(-1)**(ordinal+sum(k<i for k in rest))
                    for key,value in encode(tuple(sorted((*rest,i)))).items():
                        target[key]=s.expand(target.get(key,0)+sign*coefficient*value)
            target={key:value for key,value in target.items() if value!=0}
            assert one_site_current(source,current)==target
            checks+=1
    return {'one_object_internal_dimension':15,'two_object_exterior_dimension':105,
            'two_distinguishable_cloud_tensor_dimension':225,
            'actual_sum_of_native_currents_intertwiner_columns_checked':checks,
            'general_source_extension':'antisymmetrize whole flux-spin-position-native-cloud objects; vacuum-pad independent source slots',
            'full_exterior_CAR_encoding_is_declared_finite_architecture':True,
            'fixed_cloud_per_internal_CAR_mode_slot_is_a_current_intertwiner':False,
            'geometric_locality_of_encoded_gates_established':False,
            'native_U_on_encoded_cloud_products':'+I','native_K_on_encoded_cloud_products':'0'}


def mass_freedom():
    _,_,currents,_=representation();y=currents[-1]
    assert all(y[i,i]+y[j,j]!=0 for i in range(15) for j in range(i,15))
    blocks={q:[i for i,value in enumerate(y.diagonal()) if value==q] for q in set(y.diagonal())}
    rows=[]
    for weight,indices in sorted(blocks.items(),key=lambda row:float(row[0])):
        n=len(indices);columns=[]
        for a,b in product(range(n),repeat=2):
            unit=s.zeros(n);unit[a,b]=1
            columns.append(s.Matrix.vstack(*[(t.extract(indices,indices)*unit-unit*t.extract(indices,indices)).reshape(n*n,1)
                                              for t in currents]))
        rank=s.Matrix.hstack(*columns).rank()
        assert rank==n*n-1
        rows.append({'hypercharge':str(weight),'dimension':n,'invariant_Dirac_mass_parameters':n*n-rank})
    assert sum(row['invariant_Dirac_mass_parameters'] for row in rows)==5
    mu=s.symbols('mu0:5',real=True);mass=s.zeros(15)
    for coefficient,indices in zip(mu,blocks.values(),strict=True):
        for index in indices:mass[index,index]=coefficient
    assert all(mass*t==t*mass for t in currents)
    return {'irreducible_current_blocks':rows,'gauge_invariant_Dirac_mass_parameter_dimension':5,
            'same_unit_face_count_coin_all_internal_rest_frequencies':'2',
            'all_corresponding_internal_rest_frequency_ratios':'1',
            'constant_chiral_R15_Weyl_mass_dimension':0,
            'Dirac_R15_contains_both_chiralities_in_same_internal_representation':True,
            'face_count_coupling_controls':[{'coefficient':'1','rest_frequency':'2'},
                                            {'coefficient':'2','rest_frequency':'4'}],
            'GeV_mass_energy_formula':'hbar_in_GeV_seconds * dimensionless_rest_frequency / seconds_per_source_time_unit',
            'mass_scale_or_ratios_selected_by_charge_weights':False}
