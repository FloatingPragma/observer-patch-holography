"""Parent Haar/Slater moments retain the omitted physical Hilbert space."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

from model.particle_common_contract import contract
from model.particle_primary_variational_readout import variational
from model.relativistic_sector_audit.sector import sqrt_interval

SOURCE_INPUTS=(
    'model/particle_primary_residual_refinement/model.py',
    'model/particle_primary_residual_refinement/PROOF.md',
    'model/particle_primary_variational_readout/variational.py',
    'model/particle_primary_variational_readout/PROOF.md',
    'model/particle_primary_spectral_scale/PROOF.md',
    'model/particle_common_contract/PROOF.md',
    'model/particle_local_wilson_density/density.py',
)


def geometry(side):
    if type(side) is not int or side<2:raise ValueError('open side at least two required')
    L=side
    square_degrees=36*L**3-66*L**2+24*L
    bonds=3*L*L*(L-1)
    straight=6*L*L*(L-2)
    return {'side':L,'bonds':bonds,'plaquettes':3*L*(L-1)**2,
            'sum_squared_interior_degrees':square_degrees,
            'directed_nonbacktracking_two_step_paths':square_degrees-2*bonds,
            'directed_straight_two_step_paths':straight,
            'Haar_trace_BstarB_squared_per_native_dimension':F(square_degrees-bonds,4),
            'Haar_single_pair_second_action_norm_per_native_dimension':F(square_degrees-2*bonds+straight,4),
            'hopping_connected_cumulant_coefficient':F(square_degrees-straight,4)}


@lru_cache(maxsize=1)
def fourth_character_control():
    """Actual native character fourth moments by exact Weyl integration."""
    inherited=variational.character_control();weights=inherited['native_weights']
    order=inherited['block_order_Q6'];zero=(0,0,0,0)
    chars=[]
    for q in order:
        c=Counter(tuple(w) for w in weights if w[-1]==q)
        c.update(tuple(-x for x in w) for w in weights if w[-1]==q)
        chars.append(dict(c))
    total=Counter()
    for c in chars:total.update(c)
    weyl={zero:1}
    for root in ((2,1,0,0),(1,2,0,0),(1,-1,0,0),(0,0,2,0)):
        for sign in (-1,1):
            weyl=variational.multiply(weyl,{zero:1,tuple(sign*x for x in root):-1})
    square=variational.multiply(total,total)
    matrix=[]
    for left in chars:
        row=[]
        for right in chars:
            polynomial=variational.multiply(square,variational.multiply(left,right))
            value=F(sum(v*weyl.get(tuple(-x for x in k),0) for k,v in polynomial.items()),12)
            if value.denominator!=1:raise ArithmeticError('native Haar multiplicity is not integral')
            row.append(value.numerator)
        matrix.append(row)
    if matrix!=[list(row) for row in zip(*matrix)]:raise ArithmeticError('fourth character form is not symmetric')
    reps=[(-4,-1,0),(-3,0,1),(1,1,1),(2,-1,0),(6,0,0)]
    reps += [(-q,-color,weak) for q,color,weak in reps]
    locate={abs(q):j for j,q in enumerate(order)};independent=[[0]*5 for _ in range(5)]
    for rows in product(reps,repeat=4):
        if sum(row[0] for row in rows):continue
        fundamentals=sum(row[1]==1 for row in rows);antifundamentals=sum(row[1]==-1 for row in rows)
        colors={(0,0):1,(1,1):1,(2,2):2,(3,0):1,(0,3):1}.get((fundamentals,antifundamentals),0)
        weak={0:1,2:1,4:2}.get(sum(row[2] for row in rows),0)
        independent[locate[abs(rows[2][0])]][locate[abs(rows[3][0])]]+=colors*weak
    if independent!=matrix:raise ArithmeticError('fourth tensor singlets and Weyl integration differ')
    return {'block_order_Q6':order,'integral_S_squared_xr_xs':matrix,
            'integral_S_fourth':sum(map(sum,matrix)),
            'independent_fourfold_tensor_singlet_matrix':independent,
            'native_character_dimensions':[len([w for w in weights if w[-1]==q]) for q in order]}


def mixed_parent_moment(side,hopping_weights,loop_weights):
    """Exact <Omega,V² V_t² Omega> before link-word compression."""
    hopping_weights=tuple(hopping_weights);loop_weights=tuple(loop_weights)
    if any(type(x) not in (int,str,F) for x in hopping_weights+loop_weights):
        raise ValueError('exact rational weights required')
    th=list(map(F,hopping_weights));tl=list(map(F,loop_weights))
    if len(th)!=5 or len(tl)!=5:raise ValueError('all five native weights required')
    g=geometry(side);chars=fourth_character_control();dims=chars['native_character_dimensions']
    B=g['bonds'];P=g['plaquettes'];D=g['hopping_connected_cumulant_coefficient']
    z=15*B;zt=B*sum(d*t*t for d,t in zip(dims,th));qh=B*sum(d*t for d,t in zip(dims,th))
    zl=F(P,4410);sl=F(P,22050)*sum(t*t for t in tl);ql=F(P,22050)*sum(tl)
    hopping=z*zt+2*qh*qh-D*sum(d*t*t for d,t in zip(dims,th))
    matrix=chars['integral_S_squared_xr_xs']
    local_loop=sum(tl[r]*matrix[r][u]*tl[u] for r in range(5) for u in range(5))/210**4
    loop=P*local_loop+F(P*(P-1),210**4)*(20*sum(t*t for t in tl)+8*sum(tl)**2)
    mixed=hopping+z*sl+zt*zl+4*qh*ql+loop
    if min(hopping,loop,mixed)<0:raise ArithmeticError('a positive exact moment became negative')
    return {'geometry':g,'hopping_variance':z,'weighted_hopping_variance':zt,
            'hopping_covariance':qh,'magnetic_variance':zl,'weighted_magnetic_variance':sl,
            'magnetic_covariance':ql,'exact_hopping_mixed_fourth':hopping,
            'exact_magnetic_mixed_fourth':loop,'exact_parent_mixed_fourth':mixed,
            'full_weighted_variance':zt+sl,'full_covariance':qh+ql,
            'positive_and_negative_frequencies_or_particle_weights_computed':False}


def refined_bound(member,producer=None):
    contract.verify_member(member)
    producer=variational.bind_member(member) if producer is None else producer
    variational.verify_bound(member,producer)
    old=producer['result'];trial=old['dimensionless_trial']
    th=list(map(F,trial['rational_hopping_coefficients']))
    tl=list(map(F,trial['rational_loop_coefficients']))
    moment=mixed_parent_moment(member['source']['side'],th,tl)
    a=F(member['source']['spacing']);nl=F(trial['real_loop_packet_norm_squared'])
    q=moment['full_covariance'];square=moment['full_weighted_variance']
    if 1+square!=F(trial['trial_norm_squared_before_normalization']):
        raise ArithmeticError('the exact inherited trial normalization changed')
    correction=variational.character_control()['cubic_loop_singlet_matrix']
    projected_loop=[-sum(F(correction[i][j],420)*tl[j] for j in range(5)) for i in range(5)]
    outside=moment['exact_parent_mixed_fourth']-q*q-nl*sum(x*x for x in projected_loop)
    e=F(trial['trial_energy_shift']);dlo=-F(trial['lower_Ritz_shift_interval'][1])
    squared=((q+e)**2+(dlo+e)**2*square+outside)/(1+square)
    if outside<0 or squared<0:raise ArithmeticError('exact complement norm became negative')
    radius=sqrt_interval(squared)[1]/a
    previous=F(old['trial_residual_upper'])
    if radius>=previous:raise ArithmeticError('the exact local moments did not improve the inherited bound')
    center=F(old['trial_rayleigh_quotient'])
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'exact_parent_moments':moment,
        'inherited_trial_receipt_sha256':contract.digest(producer),
        'trial_rayleigh_quotient':center,'full_physical_ground_upper':center,
        'parent_omitted_component_squared_norm':outside,
        'dimensionless_residual_squared_upper':squared,
        'source_residual_upper':radius,
        'inherited_source_residual_upper':previous,
        'certified_radius_bound_reduction_factor':previous/radius,
        'some_full_physical_eigenvalue_interval':[center-radius,center+radius],
        'squared_residual_is_exact_on_selected_cutoff':member['source']['per_link_word_degree']>=2,
        'link_degree_one_bound_uses_contractive_final_compression':member['source']['per_link_word_degree']==1,
        'trial_subspace_invariant':False,'trial_replaces_selected_polynomial_state':False,
        'source_clock_or_reference_changed':False,
        'ground_isolation_or_particle_mass_proved':False,
        'original_primary_probe_residue_computed':False})


def bind_member(member):
    result=refined_bound(member)
    return {'result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('the actual source or full-space residual differs')


def certificate():
    return {'schema':'oph.primary_exact_local_residual.v1',
            'character_fourth_moments':fourth_character_control(),
            'members':[bind_member(contract.member(n)) for n in (0,1,2)],
            'physical_particle_masses_GeV':[]}
