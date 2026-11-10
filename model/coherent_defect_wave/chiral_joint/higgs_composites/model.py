"""Full-domain CAR polynomials with the actual five-state scalar matrix units."""
from collections import defaultdict
from functools import lru_cache
import sympy as s
from ...physical_attachment.chiral import representation
from ...chiral_higgs.model import EPS, PAULI, higgs_memory
from ..model import sequence

# Each key is (normally ordered CAR word, scalar output, scalar input).
# Creation operators precede annihilation operators; each list is increasing.
# No occupation-sector projection is used.

def clean(poly):
    return {k:s.expand(v) for k,v in poly.items() if s.expand(v)!=0}


def add(*polys):
    out=defaultdict(lambda:s.Integer(0))
    for poly in polys:
        for key,value in poly.items():out[key]+=value
    return clean(out)


def scale(poly,value):return clean({key:value*v for key,v in poly.items()})


@lru_cache(maxsize=None)
def normal_order(word):
    word=tuple(word)
    for j in range(len(word)-1):
        left,right=word[j:j+2]
        if left==right:return ()
        order=lambda atom:(not atom[1],atom[0])
        if order(left)>order(right):
            out=defaultdict(lambda:s.Integer(0))
            swapped=word[:j]+(right,left)+word[j+2:]
            for key,value in normal_order(swapped):out[key]-=value
            if not left[1] and right[1] and left[0]==right[0]:
                for key,value in normal_order(word[:j]+word[j+2:]):out[key]+=value
            return tuple(clean(out).items())
    return ((word,s.Integer(1)),)


def monomial(word=(),matrix=None):
    matrix=s.eye(5) if matrix is None else s.ImmutableSparseMatrix(matrix)
    return {(normal,r,c):value*sign for normal,sign in normal_order(tuple(word))
            for (r,c),value in s.ImmutableSparseMatrix(matrix).todok().items()}


def multiply(a,b):
    byrow=defaultdict(list)
    for (word,r,c),value in b.items():byrow[r].append((word,c,value))
    out=defaultdict(lambda:s.Integer(0))
    for (word,r,c),value in a.items():
        for other,end,weight in byrow[c]:
            for normal,sign in normal_order(word+other):out[normal,r,end]+=value*weight*sign
    return clean(out)


def adjoint(poly):
    out=defaultdict(lambda:s.Integer(0))
    for (word,r,c),value in poly.items():
        reverse=tuple((i,not create) for i,create in reversed(word))
        for normal,sign in normal_order(reverse):out[normal,c,r]+=s.conjugate(value)*sign
    return clean(out)


def commutator(a,b):return add(multiply(a,b),scale(multiply(b,a),-1))
def anticommutator(a,b):return add(multiply(a,b),multiply(b,a))


def one_body(matrix):
    return add(*(scale(monomial(((i,True),(j,False))),v)
        for (i,j),v in s.ImmutableSparseMatrix(matrix).todok().items()))


@lru_cache(maxsize=1)
def quantum_yukawa():
    """The exact three inherited pair tensors, with symbolic real coefficients."""
    fields=higgs_memory()[0]
    ys=dict(zip(('up','down','charged_lepton'),s.symbols('y_u y_d y_e',real=True),strict=True))
    parts={}
    for name,vertex in representation()[3].items():
        terms=[]
        for k,tensor in enumerate(vertex['components']):
            field=fields[k].H if vertex['conjugate_Higgs'] else fields[k]
            for (i,j),value in s.ImmutableSparseMatrix(tensor).todok().items():
                terms.append(scale(monomial(((i,False),(15+j,False)),field),ys[name]*value))
        pair=add(*terms);parts[name]=add(pair,adjoint(pair))
    return parts,ys


@lru_cache(maxsize=None)
def composite(channel,index,spin):
    """J_{i,spin}=sum_kj C^k_ij Phi_k c_{spin,j}; inherited signs retained."""
    if spin not in (0,1):raise ValueError('two Weyl components required')
    vertex=representation()[3][channel];fields=higgs_memory()[0]
    result=[]
    for k,tensor in enumerate(vertex['components']):
        phi=fields[k].H if vertex['conjugate_Higgs'] else fields[k]
        for j in range(15):
            if tensor[index,j]:result.append(scale(monomial(((15*spin+j,False),),phi),tensor[index,j]))
    return add(*result)


@lru_cache(maxsize=None)
def screened_eom(channel,index,spin):
    return add(*(scale(composite(channel,index,beta),-EPS[spin,beta]) for beta in (0,1)))


@lru_cache(maxsize=None)
def neutral_composite(spin):
    phi=higgs_memory()[0]
    return add(monomial(((15*spin+11,False),),phi[0]),
               scale(monomial(((15*spin+10,False),),phi[1]),-1))


@lru_cache(maxsize=None)
def gauge_current(current):
    matter=s.kronecker_product(s.eye(2),representation()[2][current])
    return add(one_body(matter),monomial(matrix=higgs_memory()[1][current]))


@lru_cache(maxsize=None)
def spin_current(axis):return one_body(s.kronecker_product(PAULI[axis]/2,s.eye(15)))


def column(poly,mask,scalar):
    if type(mask) is not int or not 0<=mask<2**30 or scalar not in range(5):
        raise ValueError('full native F30/scalar coordinate required')
    out=defaultdict(lambda:s.Integer(0))
    for (word,r,c),value in poly.items():
        if c!=scalar:continue
        step=sequence(mask,word)
        if step is not None:out[step[0],r]+=value*step[1]
    return clean(out)


SINGLETS={'up':(0,1,4),'down':(12,13,14),'charged_lepton':(9,)}


def tensor_certificate():
    parts,ys=quantum_yukawa();h=add(*parts.values())
    assert h==adjoint(h)
    ward=eom=screening=spin_checks=0
    for a in range(12):
        assert not commutator(gauge_current(a),h);ward+=1
    # Equality in the complete normal-ordered CAR basis proves all occupations.
    for i in range(15):
        for alpha in (0,1):
            derived=commutator(h,monomial(((15*alpha+i,True),)))
            expected=add(*(scale(screened_eom(name,i,alpha),ys[name]) for name in parts))
            assert derived==expected;eom+=1
    for channel,indices in SINGLETS.items():
        for i in indices:
            for alpha in (0,1):
                op=screened_eom(channel,i,alpha)
                for a,t in enumerate(representation()[2]):
                    expected=add(*(scale(screened_eom(channel,j,alpha),t[j,i])
                                   for j in indices if t[j,i]))
                    assert commutator(gauge_current(a),op)==expected;screening+=1
                for axis in range(3):
                    expected=add(*(scale(screened_eom(channel,i,beta),PAULI[axis][beta,alpha]/2)
                                   for beta in (0,1)))
                    assert commutator(spin_current(axis),op)==expected;spin_checks+=1
    for alpha in (0,1):
        neutral=neutral_composite(alpha)
        for a in range(12):
            assert not commutator(gauge_current(a),neutral);screening+=1
        for axis in range(3):
            expected=add(*(scale(neutral_composite(beta),-PAULI[axis][alpha,beta]/2) for beta in (0,1)))
            assert commutator(spin_current(axis),neutral)==expected;spin_checks+=1
    # Pair tensor is the SL(2,C) invariant; a finite cutoff Hamiltonian is not
    # asserted to carry boost symmetry.
    sl2=0
    for b in (s.diag(2,s.Rational(1,2)),s.Matrix([[1,s.I],[0,1]]),
              s.Matrix([[s.Rational(5,3),s.Rational(4,3)],[s.Rational(4,3),s.Rational(5,3)]])):
        assert b.det()==1 and b.T*EPS*b==EPS;sl2+=1
    return {'full_local_Fock_dimension':str(2**30),'scalar_dimension':5,
            'occupation_sectors_retained':list(range(31)),
            'full_HY_current_commutators':ward,'full_CAR_creation_EOM_identities':eom,
            'screened_current_component_identities':screening,'spin_component_identities':spin_checks,
            'SL2C_epsilon_tensor_identities':sl2,
            'native_singlet_indices':{name:list(indices) for name,indices in SINGLETS.items()},
            'operator_hypercharges':{'up':'-2/3','down':'1/3','charged_lepton':'1','neutral':'0'},
            'weak_screening':'all three weak generators vanish on the singlet composites',
            'residual_color':'up and down transform exactly as their native conjugate-singlet creation multiplets',
            'neutral_has_native_conjugate_singlet_partner':False,
            'neutral_has_Yukawa_mass_term':False,'constant_Higgs_substitution':False,
            'boost_covariance_of_finite_Hamiltonian':False}


def noncanonical_certificate():
    cases=[]
    for name,indices in SINGLETS.items():
        cases.append((name,composite(name,indices[0],0)))
    cases.append(('neutral',neutral_composite(0)))
    result=[]
    for name,op in cases:
        anti=anticommutator(op,adjoint(op))
        assert anti!=monomial()
        # Keep all five scalar alternatives and both empty/filled occupations.
        entries=[]
        for mask in (0,(1<<30)-1):
            for scalar in range(5):
                value=column(anti,mask,scalar).get((mask,scalar),s.Integer(0))
                entries.append({'fermion_occupation':'empty' if mask==0 else 'full',
                                'scalar':scalar,'diagonal_anticommutator':str(value)})
        assert len({e['diagonal_anticommutator'] for e in entries})>1
        result.append({'channel':name,'canonical_CAR':False,'complete_extreme_occupation_witnesses':entries,
                       'anticommutator_normal_ordered_terms':len(anti)})
    return result


def compiler_certificate():
    from ..periodic_family.model import yukawa_words,local_yukawa
    parts,ys=quantum_yukawa();compiled=[]
    for name,row,col,ops,coefficient in yukawa_words():
        matrix=s.zeros(5);matrix[row,col]=1
        compiled.append(scale(monomial(ops,matrix),ys[name]*coefficient))
    assert add(*compiled)==add(*parts.values())
    # Independent actual column evaluation includes every lepton occupation,
    # every scalar state, and dense occupations in all native sectors.
    active=(9,10,11,24,25,26);masks={sum(1<<active[j] for j in range(6) if m>>j&1) for m in range(64)}
    masks|={(1<<30)-1,0,sum(1<<j for j in range(30) if j%2),sum(1<<j for j in range(30) if j%3)}
    for seed in range(32):masks.add((1103515245*(seed+1)+12345)%2**30)
    h=add(*(scale(part,1) for part in parts.values()))
    substitutions={ys['up']:2,ys['down']:3,ys['charged_lepton']:5}
    numeric={key:value.subs(substitutions) for key,value in h.items()}
    count=0
    for mask in sorted(masks):
        for scalar in range(5):
            assert column(numeric,mask,scalar)==local_yukawa(mask,scalar,(s.Integer(2),s.Integer(3),s.Integer(5)))
            count+=1
    return {'exact_equality_to_periodic_full_family_Yukawa_word_compiler':True,
            'independent_full_native_column_controls':count,
            'complete_lepton_occupation_and_scalar_columns_included':320,
            'test_couplings':{'up':'2','down':'3','charged_lepton':'5'},
            'all_occupation_identity_basis':'independent normal-ordered CAR monomials times scalar matrix units'}


def composite_dynamics_certificate():
    parts,ys=quantum_yukawa();hy=add(*parts.values());number=higgs_memory()[2]
    omega=s.Symbol('omega',real=True);h=add(hy,scale(monomial(matrix=number),omega))
    rows=[]
    for name,indices in SINGLETS.items():
        for index in indices:
            op=composite(name,index,0);derivative=commutator(h,op)
            cubic={key:v for key,v in derivative.items() if len(key[0])==3}
            assert cubic
            witness=sorted(cubic,key=str)[0]
            rows.append({'channel':name,'native_index':index,'onsite_EOM_terms':len(derivative),
                         'nonzero_cubic_CAR_terms':len(cubic),
                         'cubic_witness':{'word':[[i,create] for i,create in witness[0]],
                            'scalar_row':witness[1],'scalar_col':witness[2],
                            'coefficient':str(cubic[witness])}})
    neutral=commutator(h,neutral_composite(0));assert neutral
    rows.append({'channel':'neutral','onsite_EOM_terms':len(neutral),
                 'nonzero_cubic_CAR_terms':sum(len(key[0])==3 for key in neutral),
                 'absence_of_neutral_Yukawa_partner_implies_stationarity':False})
    return {'exact_quantum_composite_EOM':rows,'linear_free_Dirac_closure':False,
            'scalar_commutators_at_cutoff_retained':True,
            'neutral_masslessness_inferred_from_operator_labels':False}


def star_selection_certificate():
    from ..field_limit.model import star_parity_checks
    inherited=star_parity_checks();checked=0
    for channel,indices in SINGLETS.items():
        for index in indices:
            for spin in (0,1):
                op=composite(channel,index,spin)
                assert op and all(len(word)%2==1 for word,_,_ in op);checked+=1
    for spin in (0,1):
        op=neutral_composite(spin)
        assert op and all(len(word)%2==1 for word,_,_ in op);checked+=1
    return {'inherited_full_domain_star_checks':inherited,
            'local_composite_star_character_checks':checked,
            'local_character':'C_x J_y C_x=(-1)^(x=y) J_y',
            'reference_hypothesis':'[rho,C_x]=0 for every x',
            'exact_separated_local_correlator':'Tr(rho J_x(t) J_y^dagger)=0 for x!=y',
            'neutral_current_singlet_evades_star_selection':False,
            'retained_link_string_dressed_correlators_excluded':False,
            'physical_dressed_neutrino_propagation_proved':False}
