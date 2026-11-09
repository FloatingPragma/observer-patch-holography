"""Full two-spin-row native gauge kernels with physical cube and current labels."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_sector_threshold import model as sectors
from model.particle_primary_charge_conjugation import model as conjugation
from model.particle_primary_channel_labels import channels
from model.cubic_particle_channels import group

OMEGA=sum(1<<j for j in range(30,60))
SOURCE_INPUTS=tuple(dict.fromkeys((
    'model/particle_primary_composite_spin/__init__.py',
    'model/particle_primary_composite_spin/__main__.py',
    'model/particle_primary_composite_spin/model.py',
    'model/particle_primary_composite_spin/test_model.py',
    'model/particle_primary_composite_spin/README.md',
    'model/particle_primary_composite_spin/PROOF.md',
    *sectors.SOURCE_INPUTS,*conjugation.SOURCE_INPUTS,*channels.SOURCE_INPUTS,
)))


def clean(A):
    # All entries lie in the fixed radical field Q(i,sqrt(2),sqrt(3)).
    return A.applyfunc(s.expand)


@lru_cache(maxsize=1)
def beta_chart():
    pairs=sectors.branch.response.beta_pairs()
    return s.Matrix.hstack(pairs[0][0],pairs[1][0],pairs[0][1],pairs[1][1])


def occupation_basis(delta):
    delta=tuple(delta)
    if len(delta)!=5 or sum(q*n for (q,_),n in zip(sectors.native_blocks(),delta)):
        raise ValueError('an actual gauge-neutral five-block increment is required')
    choices=[]
    for (_,indices),number in zip(sectors.native_blocks(),delta):
        rows=(0,1) if number>=0 else (2,3)
        modes=tuple(15*r+i for r in rows for i in indices)
        choices.append(tuple(combinations(modes,abs(number))))
    masks=[]
    for selected in product(*choices):
        mask=OMEGA
        for labels,n in zip(selected,delta):
            for i in labels:
                if n>0:mask|=1<<i
                elif n<0:mask^=1<<i
        masks.append(mask)
    return tuple(sorted(set(masks)))


def act_matrix(masks,terms):
    where={mask:i for i,mask in enumerate(masks)};A=s.zeros(len(masks))
    for column,mask in enumerate(masks):
        for target,value in sectors._act({mask:s.S.One},terms).items():
            if target not in where:raise ArithmeticError('the two-spin-row occupation domain is not invariant')
            A[where[target],column]=value
    return A


@lru_cache(maxsize=None)
def local_multiplet(delta):
    """Actual60-bit basis, orthonormal gauge kernel and Dirac-spin generators."""
    delta=tuple(delta);masks=occupation_basis(delta)
    gauges=[]
    for g in sectors.native.native_currents():
        terms=[(15*r+i,15*r+j,v) for r in range(4) for (i,j),v in g.todok().items() if v]
        gauges.append(act_matrix(masks,terms))
    kernel=s.Matrix.vstack(*gauges).nullspace()
    if not kernel:raise ArithmeticError('the local composite gauge kernel is empty')
    basis=s.Matrix.hstack(*s.GramSchmidt(kernel,orthonormal=True))
    if clean(basis.H*basis)!=s.eye(basis.cols):raise ArithmeticError('gauge-kernel basis is not orthonormal')
    if any(A*basis!=s.zeros(len(masks),basis.cols) for A in gauges):
        raise ArithmeticError('a multiplet vector is not an actual gauge singlet')
    gamma=sectors.branch.response.density.gamma_data();chart=beta_chart()
    J=[];ambient=[]
    for axis in range(3):
        spin=clean(-s.I*chart.H*gamma[(axis+1)%3]*gamma[(axis+2)%3]*chart/2)
        terms=[(15*r+i,15*t+i,spin[r,t]) for r,t in product(range(4),repeat=2)
               for i in range(15) if spin[r,t]]
        A=act_matrix(masks,terms);j=clean(basis.H*A*basis)
        if clean(A*basis-basis*j)!=s.zeros(len(masks),basis.cols):
            raise ArithmeticError('actual Dirac rotations leave the local singlet kernel')
        if j.H!=j:raise ArithmeticError('a source spin generator is not Hermitian')
        ambient.append(A);J.append(j)
    for a in range(3):
        if clean(J[a]*J[(a+1)%3]-J[(a+1)%3]*J[a]-s.I*J[(a+2)%3])!=s.zeros(basis.cols):
            raise ArithmeticError('the source spin coefficient algebra differs')
    casimir=clean(sum((j*j for j in J),s.zeros(basis.cols)))
    eigen=casimir.eigenvals();projectors={}
    for value in eigen:
        p=s.eye(basis.cols)
        for other in eigen:
            if other!=value:p=p*(casimir-other*s.eye(basis.cols))/(value-other)
        projectors[value]=clean(p)
    return {'delta':delta,'masks':masks,'gauge_matrices':tuple(gauges),
        'kernel_basis':basis,'spin_generators':tuple(J),'ambient_spin_generators':tuple(ambient),
        'spin_casimir':casimir,'casimir_multiplicities':eigen,'casimir_projectors':projectors,
        'CAR_parity':sum(map(abs,delta))%2,
        'intrinsic_inversion_parity':(-1)**sum(-n for n in delta if n<0)}


@lru_cache(maxsize=1)
def charge_mode_map():
    M=clean(beta_chart().H*sectors.branch.response.density.gamma_data()[1]*beta_chart().conjugate())
    mapping=[]
    for row in range(4):
        entries=[(j,M[row,j]) for j in range(4) if M[row,j]]
        if len(entries)!=1 or entries[0][1] not in (-1,1):
            raise ArithmeticError('the actual beta-chart charge map is not the checked signed permutation')
        for i in range(15):mapping.append((15*entries[0][0]+i,int(entries[0][1])))
    return tuple(mapping)


def charge_mask(mask):
    target=(1<<60)-1;phase=1;mapping=charge_mode_map()
    for i in reversed([j for j in range(60) if mask&(1<<j)]):
        j,c=mapping[i]
        if not target&(1<<j):raise ArithmeticError('charge map annihilated a basis vector')
        phase*=c*(-1)**((target&((1<<j)-1)).bit_count());target^=1<<j
    return target,phase


@lru_cache(maxsize=None)
def conjugate_multiplet(delta):
    delta=tuple(delta);plus=local_multiplet(delta);minus=local_multiplet(tuple(-n for n in delta))
    target={mask:i for i,mask in enumerate(minus['masks'])};U=s.zeros(len(target),len(plus['masks']))
    anchor,anchor_phase=charge_mask(OMEGA)
    if anchor!=OMEGA:raise ArithmeticError('charge conjugation does not preserve the filled-beta anchor')
    for j,mask in enumerate(plus['masks']):
        out,phase=charge_mask(mask)
        if out not in target:raise ArithmeticError('charge conjugation leaves the opposite minimal occupation domain')
        U[target[out],j]=phase*anchor_phase
    J=clean(minus['kernel_basis'].H*U*plus['kernel_basis'])
    if clean(J.H*J)!=s.eye(J.cols):raise ArithmeticError('the physical conjugate multiplets are not unitarily matched')
    for a in range(3):
        if clean(minus['spin_generators'][a]*J-J*plus['spin_generators'][a])!=s.zeros(J.rows,J.cols):
            raise ArithmeticError('conjugation fails to intertwine source spin coefficients')
    return {'ambient_map':U,'kernel_map':J,'anchor_phase':anchor_phase,
            'opposite_intrinsic_inversion':minus['intrinsic_inversion_parity']}


def defect_operators(mask):
    """Ascending quasiparticle factors, represented as (mode, creation)."""
    return tuple((i,i<30) for i in range(60) if bool(mask&(1<<i))!=(i>=30))


def apply_operators(mask,operators):
    phase=1
    for i,creation in reversed(operators):
        if bool(mask&(1<<i))==creation:return None,0
        phase*=(-1)**((mask&((1<<i)-1)).bit_count());mask^=1<<i
    return mask,phase


def composite_polynomial(delta,column):
    local=local_multiplet(tuple(delta));V=local['kernel_basis']
    if not 0<=column<V.cols:raise ValueError('multiplet column outside kernel')
    terms=[]
    for row,mask in enumerate(local['masks']):
        if V[row,column]:
            ops=defect_operators(mask);target,phase=apply_operators(OMEGA,ops)
            if target!=mask or phase not in (-1,1):raise ArithmeticError('quasiparticle sign failed')
            terms.append((ops,phase*V[row,column]))
    return tuple(terms)


@lru_cache(maxsize=None)
def exterior_rotation(delta,q):
    """Actual source spin lift acting on the full CAR gauge-kernel columns."""
    local=local_multiplet(tuple(delta));V=local['kernel_basis'];where={m:i for i,m in enumerate(local['masks'])}
    S=clean(beta_chart().H*channels.spin_lift(q)*beta_chart())
    if S[:2,2:]!=s.zeros(2) or S[2:,:2]!=s.zeros(2) or s.simplify(S[2:,2:].det())!=1:
        raise ArithmeticError('the source spin lift changes the beta anchor')
    transformed=s.zeros(V.rows,V.cols)
    for col in range(V.cols):
        for ops,coefficient in composite_polynomial(tuple(delta),col):
            choices=[]
            for mode,creation in ops:
                row,internal=divmod(mode,15)
                choices.append(tuple(((15*t+internal,creation),S[t,row] if creation else s.conjugate(S[t,row]))
                                     for t in range(4) if S[t,row]))
            for selected in product(*choices):
                out,phase=apply_operators(OMEGA,tuple(item[0] for item in selected))
                if phase:
                    transformed[where[out],col]+=coefficient*phase*s.prod(item[1] for item in selected)
    transformed=clean(transformed);D=clean(V.H*transformed)
    if clean(transformed-V*D)!=s.zeros(V.rows,V.cols) or clean(D.H*D)!=s.eye(V.cols):
        raise ArithmeticError('the actual exterior rotation does not preserve the physical multiplet')
    return D


@lru_cache(maxsize=None)
def cube_restriction(delta):
    delta=tuple(delta);local=local_multiplet(delta);matrices=[];chars=[]
    for q in group.elements():
        D=exterior_rotation(delta,q);matrices.append(D);chars.append(s.simplify(s.trace(D)))
    multiplicities={}
    for label in group.IRREPS:
        value=s.simplify(sum(channels._sym(group.characters(q)[label])*char
                             for q,char in zip(group.elements(),chars))/48)
        if not value.is_Integer or value<0:raise ArithmeticError('nonintegral actual cube character multiplicity')
        if value:multiplicities[label]=int(value)
    parity='g' if local['intrinsic_inversion_parity']==1 else 'u'
    return {'proper_binary_cube_irreps':multiplicities,
            'uniform_orbit_full_cube_irreps':{key+parity:value for key,value in multiplicities.items()},
            'characters':tuple(chars),'matrices':tuple(matrices)}


@lru_cache(maxsize=None)
def operator_controls(delta):
    """Literal quasiparticle monomials, native charges and disjoint exchange."""
    delta=tuple(delta);local=local_multiplet(delta);rows=[]
    for column in range(local['kernel_basis'].cols):
        terms=composite_polynomial(delta,column)
        for ops,_ in terms:
            counts=[]
            for _,indices in sectors.native_blocks():
                counts.append(sum((1 if creation else -1) for i,creation in ops if i%15 in indices))
            if tuple(counts)!=delta:raise ArithmeticError('a composite monomial changes the claimed native current')
        anchor=OMEGA|(OMEGA<<60)
        a=terms[0][0];b=tuple((i+60,c) for i,c in a)
        out1,p1=apply_operators(anchor,a+b);out2,p2=apply_operators(anchor,b+a)
        if out1!=out2 or p1!=(-1)**local['CAR_parity']*p2:
            raise ArithmeticError('actual disjoint-support graded exchange failed')
        rows.append({'column':column,'monomial_count':len(terms),
            'operator_norm_upper':str(s.sqrt(len(terms))),
            'anchor_success_probability':F(1,len(terms))})
    return {'native_occupation_increment':delta,'polynomial_degree':sum(map(abs,delta)),
        'gauge_invariant_composite_operator':True,'native_current_commutators_exact':True,
        'disjoint_exchange_sign':(-1)**local['CAR_parity'],
        'local_canonical_free_CAR_asserted':False,'normalized_columns':rows}


def multiplet_record(delta):
    local=local_multiplet(tuple(delta));cube=cube_restriction(tuple(delta));C=conjugate_multiplet(tuple(delta))
    return {'delta':delta,'local_occupation_dimension':len(local['masks']),
        'physical_gauge_kernel_dimension':local['kernel_basis'].cols,'actual_gauge_generators':12,
        'spin_coefficient_Casimir_multiplicities':{str(k):int(v) for k,v in local['casimir_multiplicities'].items()},
        'proper_binary_cube_irreps':cube['proper_binary_cube_irreps'],
        'uniform_orbit_full_cube_irreps':cube['uniform_orbit_full_cube_irreps'],
        'CAR_parity':local['CAR_parity'],'intrinsic_inversion_parity':local['intrinsic_inversion_parity'],
        'literal_charge_conjugation_unitary_intertwiner':True,'local_C_anchor_phase':C['anchor_phase'],
        'opposite_delta_intrinsic_inversion':C['opposite_intrinsic_inversion'],
        'composite_operator_controls':operator_controls(tuple(delta)),
        'local_kernel_H_invariance_asserted':False,'continuous_rotation_or_Lorentz_spin_asserted':False}


def bind_member(member):
    sectors.branch._member_zero(member)
    bound=sectors.bind_member(member);sectors.verify_bound(member,bound)
    conjugate_bound=conjugation.bind_member(member);conjugation.verify_bound(member,conjugate_bound)
    rows=[];ground_lower=F(bound['result']['global_ground_lower']);ground_upper=F(bound['result']['global_ground_upper'])
    for family,(_,delta,_) in zip(bound['result']['sector_results'],sectors.SECTORS):
        signed=[]
        for sector in family['signed_sectors']:
            actual=tuple(int(n) for n in sector['delta']);lower=F(sector['physical_sector_energy_lower'])
            mean=F(-720+3*sum(map(abs,actual)))
            record=sectors.prepared_record(member,{'sector_ground_energy_upper':mean},lower,ground_lower,ground_upper)
            record['source_post_intervention_state']='arbitrary normalized uniform-orbit local gauge-singlet spin multiplet'
            record['complete_cube_projective_readout_commutes_with_H']=True
            record['all_cube_irreps_and_energy_records_retained']=True
            signed.append({'multiplet':multiplet_record(actual),'bare_source_energy_mean':mean,
                'physical_full_sector_energy_lower':lower,
                'each_listed_cube_sector_bottom_excitation_interval':[lower-ground_upper,mean-ground_lower],
                'prepared_joint_cube_energy_readout':record})
        rows.append({'label':family['label'],'signed_multiplets':signed})
    result=contract.serial({'member_id':member['member_id'],**member['ids'],'native_charge_order':[-4,-3,1,2,6],
        'source_sector_receipt_sha256':contract.digest(bound),'charge_conjugation_receipt_sha256':contract.digest(conjugate_bound),
        'multiplet_families':rows,'original_state_retained_with_reference_in_idle_memory':True,
        'source_preparation_changes_selected_marginal':True,'same_source_action_reference_clock':True,
        'named_standard_model_species_identified':False,'physical_mass_table_GeV':[]})
    return {'schema':'oph.primary_composite_spin.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('physical composite spin or preparation result changed')


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
    actual={'schema':'oph.primary_composite_spin_certificate.v1','selected_member_bound':bind_member(contract.member(0)),
        'literal_source_rotation_control':channels.rotation_control(),
        'paired_memory_control':sectors.paired_memory_control(),'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('source inputs changed during composite-spin replay')
    return actual
