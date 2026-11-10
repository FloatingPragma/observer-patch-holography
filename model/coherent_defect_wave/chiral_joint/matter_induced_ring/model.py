"""Exact physical projection, complete matter elimination and loop ancestry."""
from collections import defaultdict
from functools import lru_cache
from itertools import product
import sympy as s
from ..model import active_terms,sequence
from ..two_cell.basis import current_data
from ..two_cell.symmetries import periodic_column,periodic_current_column
from ..periodic_spectrum.basis import coordinate_rows,character_multiplicities
from ..periodic_spectrum.spectrum import resolved_block
from ..periodic_boson_records.model import operators

BLOCK=(-1,-1,1,0)
MATTER_MASK=64
PARTS=('diagonal','hop01','hop10','Y0','Y1')
P_NAMES=('blank','scalar_link_10','singlet_loop','scalar_link_01','doublet_loop')
NORM_BOUND=s.Integer(30)


@lru_cache(maxsize=1)
def physical_split():
    keys,embedding,g,k,h=resolved_block(*BLOCK)
    supports=defaultdict(set)
    for (i,j),value in embedding.todok().items():supports[j].add(keys[i][0])
    p=[];q=[]
    for j in range(h.rows):
        assert MATTER_MASK not in supports[j] or supports[j]=={MATTER_MASK}
        (p if supports[j]=={MATTER_MASK} else q).append(j)
    assert (len(keys),len(p),len(q))==(345,5,67)
    assert h.rows==character_multiplicities(-1)[-1,1]==len(p)+len(q)
    assert g.extract(p,q)==s.zeros(len(p),len(q))
    gp=g.extract(p,p);gq=g.extract(q,q)
    normalization=s.diag(1,s.Rational(1,2),1,s.Rational(1,2),s.Rational(1,2))
    assert normalization.T*gp*normalization==s.eye(5)
    a=normalization.inv()*h.extract(p,p)*normalization
    assert a==s.diag(0,3,2,3,2)
    d=h.extract(q,q);b=h.extract(q,p)*normalization;c=normalization.T*gp*h.extract(p,q)
    assert c==b.T*gq and d.T*gq==gq*d
    ep=embedding[:,p]*normalization
    assert ep.H*ep==s.eye(5)
    seed=s.zeros(len(keys),1);seed[keys.index((64,0,0,0,0))]=1
    assert ep[:,0]==seed
    for name,index in (('singlet_loop',2),('doublet_loop',4)):
        op=operators(name)[0];raw=s.zeros(len(keys),1)
        for (row,col),value in op.todok().items():
            if col==0:
                left,right=divmod(row,6);raw[keys.index((64,left,right,0,0))]=value
        assert ep[:,index]==raw
    return {'keys':keys,'embedding':embedding,'g':g,'k':k,'h':h,
            'p':tuple(p),'q':tuple(q),'gp':gp,'gq':gq,'normalization':normalization,
            'a':a,'d':d,'b':b,'c':c,'ep':ep}


@lru_cache(maxsize=1)
def exact_parts():
    block=physical_split();keys=block['keys'];positions={key:i for i,key in enumerate(keys)}
    phases=[s.I**((key[0]&63).bit_count()) for key in keys]
    rows={name:{} for name in PARTS}
    for j,key in enumerate(keys):
        for target,value in periodic_column(key).items():
            assert target in positions
            differences=[n for n in range(1,5) if key[n]!=target[n]]
            if target==key:name='diagonal';assert not differences
            else:
                assert len(differences)==1
                name={1:'hop01',2:'hop10',3:'Y0',4:'Y1'}[differences[0]]
            i=positions[target]
            rows[name][i,j]=s.expand(s.conjugate(phases[i])*value*phases[j])
    rational_embedding=s.ImmutableSparseMatrix(len(keys),72,{
        (i,j):s.expand(s.conjugate(phases[i])*value)
        for (i,j),value in block['embedding'].todok().items()})
    parts={}
    for name,entries in rows.items():
        raw=s.ImmutableSparseMatrix(len(keys),len(keys),entries)
        action=raw*rational_embedding
        restricted=s.ImmutableSparseMatrix(action.extract(coordinate_rows(-1,-1,1),range(72)))
        assert rational_embedding*restricted==action
        assert restricted.T*block['g']==block['g']*restricted
        parts[name]=restricted
    assert sum(parts.values(),s.zeros(72))==block['h']
    return parts


def projection_current_checks():
    """The matter projector commutes with all currents on the full CAR domain."""
    matter,_,_=current_data();count=0
    for local_mask in (0,1):
        p=s.zeros(64);p[local_mask,local_mask]=1
        for current in matter:
            entries=defaultdict(lambda:s.Integer(0))
            for mask in range(64):
                for (i,j),value in current.todok().items():
                    for spin in range(2):
                        target=sequence(mask,((3*spin+i,True),(3*spin+j,False)))
                        if target:entries[target[0],mask]+=value*target[1]
            q=s.ImmutableSparseMatrix(64,64,dict(entries))
            assert p*q==q*p;count+=1
    return count


def gauss_embedding_checks():
    block=physical_split();count=0
    for site,current in product(range(2),range(12)):
        images=defaultdict(lambda:s.Integer(0))
        for (i,j),coefficient in block['embedding'].todok().items():
            for target,value in periodic_current_column(block['keys'][i],site,current,(0,-1)).items():
                images[target,j]+=coefficient*value
        assert all(s.expand(value)==0 for value in images.values());count+=1
    return count


def norm_certificate():
    _,cs,fields,labels=active_terms()
    yukawa_terms=2*sum(len(c.todok()) for c in cs)
    assert yukawa_terms==8 and len(labels)==6
    assert all(abs(value)==1 for c in cs for value in c.todok().values())
    for field in fields:
        a=field.H*field
        assert a*a==a and a==a.H
    hopping_terms=2*2*2*(len(labels)-1)
    assert 2*yukawa_terms+4+s.Rational(hopping_terms,4)==NORM_BOUND
    return {'full_periodic_carrier_dimension':2**12*6**2*5**2,
            'local_unit_norm_Yukawa_terms_per_cell':yukawa_terms,
            'unit_norm_scalar_and_link_projectors':4,
            'quarter_norm_oriented_hopping_terms':hopping_terms,
            'full_Hamiltonian_norm_upper_bound':str(NORM_BOUND),
            'all_occupation_scalar_link_alternatives_covered':True}


def split_certificate():
    block=physical_split();parts=exact_parts();p=block['p'];ep=block['ep']
    for name in PARTS[1:]:assert parts[name].extract(p,p)==s.zeros(5)
    initial=s.zeros(72,1);initial[p[0]]=1
    assert parts['diagonal']*initial==s.zeros(72,1)
    states=[]
    for j,name in enumerate(P_NAMES):
        states.append({'name':name,'energy':str(block['a'][j,j]),
            'normalized_raw_coordinates':[{'key':list(block['keys'][i]),'amplitude':str(value)}
                for (i,c),value in sorted(ep.todok().items()) if c==j]})
    return {'boundary_hypercharge_compensator':[0,-1],'lepton_number':-1,'two_Sz':1,
            'complete_diagonal_Gauss_weight_dimension':345,'complete_physical_dimension':72,
            'retained_matter_mask':MATTER_MASK,'retained_dimension':5,'eliminated_dimension':67,
            'retained_coordinate_indices':list(p),'eliminated_coordinate_indices':list(block['q']),
            'retained_orthonormal_states':states,
            'full_CAR_projector_current_commutators':projection_current_checks(),
            'complete_embedding_endpoint_Gauss_identities':gauss_embedding_checks(),
            'orthogonal_P_Q_Gram_blocks':True,'complete_H_and_each_part_invariance':True,
            'PHP_contains_only_diagonal_scalar_and_link_energies':True,
            'full_electric_scalar_only_H_leaves_initial_state_fixed':True,
            'P_plus_Q_is_the_full_conserved_physical_block':True,
            'quantum_Yukawa_and_all_eliminated_matter_states_retained':True}


def irreducible_words(length,target):
    """Chronological actual-H terms, with Q after every intermediate action."""
    block=physical_split();parts=exact_parts();p=block['p'];g=block['g']
    seed=s.zeros(72,1);seed[p[0]]=1
    bra=s.zeros(1,72);bra[0,p[target]]=block['normalization'][target,target]
    bra=bra*g
    states={():seed}
    for step in range(length):
        following={}
        for word,vector in states.items():
            for name,matrix in parts.items():
                image=s.MutableDenseMatrix(matrix*vector)
                if step<length-1:
                    for i in p:image[i]=0
                if image!=s.zeros(72,1):following[word+(name,)]=image
        states=following
    result={word:s.factor((bra*vector)[0]) for word,vector in states.items() if (bra*vector)[0]}
    expected=(block['c']*block['d']**(length-2)*block['b'])[target,0]
    assert sum(result.values(),s.Integer(0))==expected
    return result


def ancestry_certificate():
    block=physical_split();b,c,d=block['b'],block['c'],block['d'];moments=[]
    for n in range(6):
        kernel=c*(d**n)*b
        assert kernel==kernel.T
        moments.append({'Q_internal_powers':n,'total_H_factors':n+2,
            'blank_to_singlet_loop':str(kernel[2,0]),'blank_to_doublet_loop':str(kernel[4,0]),
            'blank_return_self_energy':str(kernel[0,0])})
    assert [row['blank_to_singlet_loop'] for row in moments[:2]]==['1/16','1/16']
    assert [row['blank_to_doublet_loop'] for row in moments[:3]]==['0','0','3/8']
    words=[]
    for length,target in ((2,2),(4,4)):
        result=irreducible_words(length,target)
        for word,value in result.items():
            assert word.count('hop01')==word.count('hop10')==1
            if length==4:assert word.count('Y0')+word.count('Y1')==2
            words.append({'target':P_NAMES[target],'chronological_terms':list(word),'coefficient':str(value)})
    # The two-edge term is absent if either actual edge part is removed.
    parts=exact_parts();p,q=block['p'],block['q'];normal=block['normalization']
    for omitted in ('hop01','hop10'):
        reduced=sum((m for name,m in parts.items() if name!=omitted),s.zeros(72))
        v=reduced.extract(q,p)*normal
        w=normal.T*block['gp']*reduced.extract(p,q)
        assert (w*v)[2,0]==0
    return {'complete_Q_irreducible_moments':moments,'nonzero_minimal_term_histories':words,
            'same_edge_return_terms_retained_in_full_self_energy':True,
            'singlet_ring_needs_both_distinct_oriented_edges':True,
            'doublet_minimal_ring_has_two_quantum_Yukawa_factors':True,
            'irreducibility_means_no_intermediate_P_return':True,
            'two_edge_noncontractible_holonomy_is_a_plaquette':False}


@lru_cache(maxsize=None)
def eliminated_kernel(z):
    z=s.sympify(z)
    if z.is_Rational is not True or z>=-NORM_BOUND:
        raise ValueError('exact rational z < -30 required by the declared resolvent certificate')
    block=physical_split();p,q=block['p'],block['q']
    gq=block['gq'];d=block['d'];b=block['b'];a=block['a']
    denominator=s.ImmutableSparseMatrix(gq*(d-z*s.eye(67)))
    coupling=gq*b
    assert denominator==denominator.T
    solved=s.ImmutableSparseMatrix(denominator.inv(method='DM')*coupling)
    assert denominator*solved==coupling
    effective=s.ImmutableSparseMatrix(a-coupling.T*solved)
    assert effective==effective.T
    lift=s.zeros(72,5)
    for j,index in enumerate(p):lift[index,:]=block['normalization'][j,:]
    for j,index in enumerate(q):lift[index,:]=-solved[j,:]
    residual=(block['h']-z*s.eye(72))*lift
    target=s.zeros(72,5)
    for j,index in enumerate(p):target[index,:]=(block['normalization']*(effective-z*s.eye(5)))[j,:]
    assert residual==target
    return effective


def resolvent_certificate():
    evaluated=[]
    for z in (-32,-64):
        kernel=eliminated_kernel(z)
        assert kernel[2,0]<0 and kernel[4,0]<0
        evaluated.append({'z':str(z),'certified_distance_to_Q_spectrum_lower_bound':str(-z-NORM_BOUND),
            'blank_to_singlet_loop':str(kernel[2,0]),'blank_to_doublet_loop':str(kernel[4,0]),
            'complete_effective_five_by_five_matrix':[[str(v) for v in kernel.row(i)] for i in range(5)],
            'exact_full_72_state_Feshbach_reconstruction_identity':True})
    assert eliminated_kernel(-32)[2,0]!=eliminated_kernel(-64)[2,0]
    assert eliminated_kernel(-32)[4,0]!=eliminated_kernel(-64)[4,0]
    return {'definition':'H_eff(z)=A-C(D-zI)^(-1)B; F(z)=H_eff(z)-zI',
            'certified_real_domain':'z < -30','evaluations':evaluated,
            'energy_dependent_loop_coefficients_verified':True,
            'all_P_and_Q_states_in_exact_inverse':True,'electric_only_off_diagonal_loop_entries':['0','0'],
            'large_negative_z_singlet_leading_coefficient':'1/(16 z)',
            'large_negative_z_doublet_leading_coefficient':'3/(8 z^3)',
            'autonomous_energy_independent_pure_gauge_generator_derived':False}
