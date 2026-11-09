"""Complete G1/H creation probes and exact full-CAR energy-moment bounds."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import lcm

import numpy as np
import sympy as s
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components

from model.particle_common_contract import contract, provenance
from model.particle_primary_composite_spin import model as spin
from model.particle_primary_selected_odd_response import model as odd
from model.particle_primary_charged_readout import model as meters
from model.particle_primary_pair_channels import model as pair_channels
from model.particle_primary_composite_transport import model as transport

DELTA=(0,2,0,0,1)
SOURCE_INPUTS=tuple(dict.fromkeys((
 'model/particle_primary_selected_spin_response/__init__.py',
 'model/particle_primary_selected_spin_response/__main__.py',
 'model/particle_primary_selected_spin_response/model.py',
 'model/particle_primary_selected_spin_response/test_model.py',
 'model/particle_primary_selected_spin_response/README.md',
 'model/particle_primary_selected_spin_response/PROOF.md',
 *spin.SOURCE_INPUTS,*odd.SOURCE_INPUTS,*meters.SOURCE_INPUTS,*pair_channels.SOURCE_INPUTS,*transport.SOURCE_INPUTS,
)))


def clean(A):return A.applyfunc(s.expand)


@lru_cache(maxsize=1)
def local_controls():
    local=spin.local_multiplet(DELTA);blocks=dict(spin.sectors.native_blocks())
    active=tuple(sorted(15*r+i for r in (0,1) for i in blocks[-3]+blocks[6]));where={j:i for i,j in enumerate(active)}
    original=[[(tuple(where[i] for i,creation in ops),value)
               for ops,value in spin.composite_polynomial(DELTA,col)] for col in range(6)]
    number=s.diag(*[mask.bit_count() for mask in range(64)])
    families={};total=s.zeros(64)
    for label,casimir in (('G1',s.Rational(3,4)),('H',s.Rational(15,4))):
        P=local['casimir_projectors'][casimir]
        W=s.Matrix.hstack(*s.GramSchmidt(P.columnspace(),orthonormal=True));polynomials=[];matrices=[]
        for col in range(W.cols):
            coefficients={}
            for k in range(6):
                for word,value in original[k]:
                    coefficients[word]=s.expand(coefficients.get(word,0)+value*W[k,col])
            coefficients={word:value for word,value in coefficients.items() if value}
            polynomials.append(coefficients)
            matrices.append(clean(sum((value*odd.creation_word(word,6) for word,value in coefficients.items()),s.zeros(64))))
        effect=clean(sum((A.H*A for A in matrices),s.zeros(64)))
        lower=clean(effect-W.cols*s.eye(64)+s.Rational(W.cols,2)*number)
        if any(value<0 for value in lower.eigenvals()):raise ArithmeticError('spin-multiplet occupation inequality failed')
        if any(value<0 for value in (W.cols*s.eye(64)-effect).eigenvals()):raise ArithmeticError('spin success cap failed')
        if any(any(value<0 for value in (s.eye(64)-A.H*A).eigenvals()) for A in matrices):
            raise ArithmeticError('a normalized weak composite is not a contraction')
        total+=effect
        families[label]={'dimension':W.cols,'coefficient_basis':W,'polynomials':tuple(polynomials),
            'operators':tuple(matrices),'effect':effect,'effect_spectrum':effect.eigenvals(),
            'deficit_spectrum':lower.eigenvals()}
    if any(value<0 for value in (6*s.eye(64)-total).eigenvals()):raise ArithmeticError('six-row instrument is not complete')
    return {'active_source_modes':active,'positive_number':number,'families':families,'total_success_effect':total/6}


def _integer_polynomial(polynomial):
    first=next(iter(polynomial.values()));ratios={w:s.simplify(v/first) for w,v in polynomial.items()}
    if not all(value.is_Rational for value in ratios.values()):raise ArithmeticError('common radical coefficient absent')
    den=lcm(*(int(value.q) for value in ratios.values()))
    square=s.simplify((first/den)**2)
    if not square.is_Rational or square<=0:raise ArithmeticError('positive rational polynomial weight required')
    return F(square),{w:int(value*den) for w,value in ratios.items()}


@lru_cache(maxsize=64)
def _creation_sparse(word,modes):
    rows=[];columns=[];values=[]
    for mask in range(1<<modes):
        target,phase=spin.apply_operators(mask,tuple((i,True) for i in word))
        if phase:rows.append(target);columns.append(mask);values.append(phase)
    return coo_matrix((np.array(values,dtype=np.int64),(rows,columns)),shape=(1<<modes,1<<modes)).tocsr()


@lru_cache(maxsize=2)
def residual_gram(label):
    family=local_controls()['families'][label];rows=[]
    for polynomial in family['polynomials']:
        square,coefficients=_integer_polynomial(polynomial);R=csr_matrix((4096,4096),dtype=np.int64)
        for word,coefficient in coefficients.items():
            for pos in range(3):
                modified=list(word);modified[pos]+=6
                R+=coefficient*_creation_sparse(tuple(modified),12)
        rows.append((square,R))
    denominator=lcm(*(square.denominator for square,_ in rows))
    G=sum((int(denominator*square)*(R.T@R) for square,R in rows),csr_matrix((4096,4096),dtype=np.int64))
    G.eliminate_zeros()
    if (G-G.T).nnz:raise ArithmeticError('residual Gram is not exactly symmetric')
    return denominator,G


@lru_cache(maxsize=2)
def residual_control(label):
    den,G=residual_gram(label);upper=F(19,3) if label=='G1' else F(12)
    count,labels=connected_components(G,directed=False);minimum=None;zero=0;largest=0
    for component in range(count):
        indices=np.flatnonzero(labels==component);largest=max(largest,len(indices))
        A=s.Rational((upper*den).numerator,(upper*den).denominator)*s.eye(len(indices))-s.Matrix(G[indices,:][:,indices].toarray())
        if A==s.zeros(len(indices)):zero+=1;continue
        L,D=A.LDLdecomposition(hermitian=False)
        if L*D*L.T!=A:raise ArithmeticError('exact rational residual factorization failed')
        pivots=list(D.diagonal())
        if any(value<=0 for value in pivots):raise ArithmeticError('the residual norm certificate is not positive')
        smallest=min(pivots);minimum=smallest if minimum is None or smallest<minimum else minimum
    counter=None
    if label=='G1':
        masks=(224,336,392,644,770,1092,1281,2114,2177)
        vector=s.Matrix((-55,-32,32,-28,8,28,-8,40,-40))
        matrix=s.Matrix(G[list(masks),:][:,list(masks)].toarray())
        numerator=int((vector.T*matrix*vector)[0]);normalizer=den*int((vector.T*vector)[0])
        if numerator<=6*normalizer:raise ArithmeticError('the full-CAR 3d countercontrol failed')
        counter={'masks':masks,'integer_vector':[int(v) for v in vector],'Rayleigh_numerator':numerator,
                 'Rayleigh_denominator':normalizer,'strictly_above_six':True,
                 'domain':'uncompressed twelve-mode CAR residual; not a physical-sector eigenvalue'}
    return {'label':label,'normalized_hopping_image_modes':6,'full_CAR_dimension':4096,
        'integer_Gram_denominator':den,'connected_blocks':count,'largest_exact_block':largest,
        'zero_blocks_at_upper_bound':zero,'smallest_positive_LDL_pivot':F(minimum),
        'squared_residual_norm_upper':upper,'exact_rational_factorization_checked':True,
        'three_times_dimension_shortcut_countercontrol':counter}


@lru_cache(maxsize=1)
def covariance_controls():
    local=spin.local_multiplet(DELTA);control=local_controls();cube=spin.cube_restriction(DELTA)
    rows={}
    for label,family in control['families'].items():
        W=family['coefficient_basis'];representations=[]
        for D in cube['matrices']:
            projected=clean(W.H*D*W)
            if clean(D*W-W*projected)!=s.zeros(6,W.cols):raise ArithmeticError('probe irrep not source-covariant')
            representations.append(projected)
        character=[s.trace(D) for D in representations]
        squared=s.simplify(sum(s.conjugate(x)*x for x in character)/48)
        if squared!=1:raise ArithmeticError('the spin-resolved probe is not irreducible')
        rows[label]={'dimension':W.cols,'actual_binary_rotations':48,'character_inner_product':1}
    minus=spin.conjugate_multiplet(DELTA)
    if clean(minus['kernel_map'].H*minus['kernel_map'])!=s.eye(6):raise ArithmeticError('source C map differs')
    return {'probe_irreps':rows,'complete_probe_rows_per_site':6,
        'source_native_increment':DELTA,'source_conjugate_increment':tuple(-n for n in DELTA),
        'irreducible_probe_label_is_final_total_state_label':False,
        'component_energy_laws_equal_after_uniform_address_average':True,
        'initial_selected_state_is_assumed_a_cube_singlet':False,
        'C_conjugate_complete_record_laws_equal':True}


@lru_cache(maxsize=1)
def current_meter_control():
    control=local_controls();active=control['active_source_modes'];blocks=spin.sectors.native_blocks()
    numbers=[s.diag(*[sum(bool(mask&(1<<j)) for j,i in enumerate(active) if i%15 in indices)
                      for mask in range(64)]) for _,indices in blocks]
    for family in control['families'].values():
        for A in family['operators']:
            for N,delta in zip(numbers,DELTA):
                if N*A-A*N!=delta*A:raise ArithmeticError('native current increment differs')
    decoded={tuple((sign*d)%3 for d in DELTA):tuple(sign*d for d in DELTA) for sign in (0,1,-1)}
    if len(decoded)!=3:raise ArithmeticError('native transfer alphabet aliases')
    for initial in product(range(3),repeat=5):
        for sign in (0,1,-1):
            after=tuple((n+sign*d)%3 for n,d in zip(initial,DELTA))
            record=tuple((a-b)%3 for a,b in zip(after,initial))
            if decoded[record]!=tuple(sign*d for d in DELTA):raise ArithmeticError('coherent cancellation depends on initial occupation')
    return {'native_occupation_order':[-4,-3,1,2,6],'radix':3,'meter_count':5,
        'complete_pointer_dimension':243,'transfer_alphabet':[tuple(-d for d in DELTA),(0,)*5,DELTA],
        'native_current_commutators_checked':True,'initial_absolute_numbers_measured':False,
        'two_coupling_identity_preserves_full_source_reference_histories':True,
        'extra_transfers_outside_alphabet_could_alias':True}


def spatial_analyzer_control():
    cells=tuple(product(range(2),repeat=3))
    W=s.Matrix([[(-1)**sum(a*b for a,b in zip(k,x)) for x in cells] for k in cells])/s.sqrt(8)
    if clean(W.H*W)!=s.eye(8):raise ArithmeticError('complete spatial pointer rotation is not unitary')
    return {'member_zero_spatial_matrix':'eight-row Walsh unitary on the coherent success address',
        'coherent_address_is_not_recorded_before_pointer_rotation':True,
        'all_spatial_rows_retained':True,'aggregate_each_internal_spin_channel_weight_unchanged':True,
        'a_fixed_spatial_row_positive_weight_proved':False,'spatial_label_is_conserved_momentum':False,
        'prepared_anchor_transport_moment_transferred_to_selected_spectrum':False}


def bind_member(member):
    contract.verify_member(member);inherited=odd.bind_member(member);odd.verify_bound(member,inherited)
    data=inherited['result'];a=F(member['source']['spacing']);V=member['source']['side']**3
    defect=F(data['weak_plus_e_defect_count_upper']);vacancy=1-defect/(4*V)
    if vacancy<F(59,400):raise ArithmeticError('selected vacancy floor changed')
    degree=F(2*member['source']['interior_hopping_bonds'],V)
    plan=odd.spectroscopy.clock_plan(member,F(1,8)/a,'1/16384','1/16384')
    clock_loss=2*F(plan['clock_plan']['complete_clock_tail_upper'])+F(plan['two_complete_clock_error_upper'])
    if clock_loss>F(1,4096):raise ArithmeticError('complete spin-clock error exceeded its budget')
    rows=[]
    for label,family in local_controls()['families'].items():
        dimension=family['dimension'];C=residual_control(label)['squared_residual_norm_upper']
        probability=F(dimension,6)*vacancy;moment=C*degree/(12*a*a)
        window=(F(1,2)/a,F(35,2)/a);weight=probability-moment/(F(17,2)/a)**2
        simple=F(1,200) if label=='G1' else F(3,200)
        if weight-clock_loss<=simple:raise ArithmeticError('positive original spin-resolved clock event failed')
        sharp=None
        if member['source']['level']==0:
            narrow=probability-moment/36
            simple_narrow=F(121,512) if label=='G1' else F(61,128)
            if narrow-clock_loss<=simple_narrow:raise ArithmeticError('member-zero spin-resolved event failed')
            sharp={'original_transfer_interval':[3,15],'original_weight_lower':narrow,
                'finite_reported_transfer_interval':[F(11,4),F(61,4)],
                'finite_reported_original_weight_lower':narrow-clock_loss,
                'simple_reported_original_weight_lower':simple_narrow,
                'final_excitation':odd.readout.excitation_window((F(3),F(15)),narrow,F(member['state']['target_ground_excess']),F(1))}
        rows.append({'probe_irrep':label,'dimension':dimension,'source_current_increment':DELTA,
            'original_success_weight_lower':probability,'second_moment_center':9/a,
            'centered_second_moment_upper':moment,'original_transfer_interval':window,
            'original_transfer_weight_lower':weight,'finite_record_interval':[F(1,4)/a,F(71,4)/a],
            'finite_record_original_weight_lower':weight-clock_loss,'simple_reported_weight_lower':simple,
            'each_component_reported_weight_lower':(weight-clock_loss)/dimension,
            'weights_are_original_and_not_success_normalized':True,'member_zero_sharp_record':sharp})
    result=contract.serial({'member_id':member['member_id'],**member['ids'],
        'selected_odd_receipt_sha256':contract.digest(inherited),'positive_beta_defect_vacancy_lower':vacancy,
        'source_site_count':V,'average_graph_degree':degree,'complete_local_outcomes_per_address':7,
        'complete_address_and_success_failure_outcomes':7*V,'complete_clock_plan':plan,
        'complete_two_clock_probability_error_upper':clock_loss,'probe_channels':rows,
        'coherent_native_current_meter':current_meter_control(),'covariance':covariance_controls(),
        'spatial_analyzer':spatial_analyzer_control(),'source_input_is_original_selected_density':True,
        'source_preparation_replaced_by_anchor_or_trial':False,'all_failures_memory_and_reference_retained':True,
        'compensating_memory_native_increment':tuple(-d for d in DELTA),
        'total_source_memory_native_numbers_and_CAR_parity_preserved':True,
        'fixed_positive_energy_survival_under_refinement_proved':False,
        'Lorentz_spin_or_named_particle_identified':False,'physical_mass_table_GeV':[]})
    return {'schema':'oph.primary_selected_spin_response.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('selected spin response or source binding changed')


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
    local=local_controls();summary={label:{'dimension':row['dimension'],
        'effect_spectrum':{str(k):int(v) for k,v in row['effect_spectrum'].items()},
        'deficit_spectrum':{str(k):int(v) for k,v in row['deficit_spectrum'].items()},
        'polynomials':[[{'modes':word,'coefficient':str(c)} for word,c in p.items()] for p in row['polynomials']]}
        for label,row in local['families'].items()}
    value={'schema':'oph.primary_selected_spin_response_certificate.v1',
        'exact_local_controls':summary,'exact_residual_controls':contract.serial([residual_control(k) for k in ('G1','H')]),
        'members':[bind_member(contract.member(n)) for n in (0,1,2,8)],'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('source changed during selected spin replay')
    return contract.serial(value)
