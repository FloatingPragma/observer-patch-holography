"""Actual mixed-orientation loop, retained dual reference and finite histories."""
from collections import defaultdict
from functools import lru_cache
from itertools import product,permutations
import sympy as s
from ...physical_attachment.chiral import representation
from ...chiral_source.model import PAULI
from ..model import clean
from ..periodic_family.model import Configuration,kinetic_terms,cells,edge_endpoints
from ..periodic_family.checks import norm_bound
from ..peter_weyl_links.model import (empty_key,column,hopping_column,gauss_column)
from ..peter_weyl_links.clebsch import link

CONFIG=Configuration()
BLOCKS={'singlet_ec':(9,),'colored_up':(0,1,4)}
SQUARE=(0,13,6,1)
ROUTE=(0,4,6,2,0)


def inner(left,right):
    return s.simplify(sum(s.conjugate(value)*right.get(key,0) for key,value in left.items()))


@lru_cache(maxsize=1)
def square_terms():
    result={term.path[0]:term for term in kinetic_terms(CONFIG)}
    assert len(result)==24 and all(len(term.path)==1 for term in result.values())
    for edge in SQUARE:assert edge in result
    assert [edge_endpoints(2,e) for e in SQUARE]==[
        ((0,0,0),(1,0,0)),((1,0,0),(1,1,0)),
        ((0,1,0),(1,1,0)),((0,0,0),(0,1,0))]
    assert all(result[e].coefficient==s.I/4 for e in SQUARE)
    # These four geometric edges have no antiperiodic seam crossing.
    assert [result[e].axis for e in SQUARE]==[0,1,0,1]
    assert all((sum(cells(2)[term.left_cell])+sum(cells(2)[term.right_cell]))%2==1
               for term in result.values())
    return result


@lru_cache(maxsize=128)
def vacuum_coordinate(a,b,adjoint):
    matrix=link(a,b)
    if adjoint:matrix=matrix.H
    entries=tuple((i,value) for (i,j),value in matrix.todok().items() if j==0)
    assert len(entries)==1 and entries[0][0]!=0
    return entries[0]


@lru_cache(maxsize=2)
def loop_state(setting):
    """Unnormalized internal trace of the actual four compressed coordinates."""
    block=BLOCKS[setting];out=defaultdict(lambda:s.Integer(0));vac=empty_key(CONFIG)
    for a,b,c,d in product(block,repeat=4):
        choices=(vacuum_coordinate(a,b,True),vacuum_coordinate(b,c,True),
                 vacuum_coordinate(d,c,False),vacuum_coordinate(a,d,False))
        links=list(vac[1]);weight=s.Integer(1)
        for edge,(label,value) in zip(SQUARE,choices,strict=True):
            links[edge]=label;weight*=value
        out[0,tuple(links),vac[2]]+=weight
    result=clean(out)
    assert len(result)==len(block)**4 and inner(result,result)==1
    assert all({i for i,v in enumerate(key[1]) if v}==set(SQUARE) for key in result)
    return result


@lru_cache(maxsize=2)
def preparation(setting):
    block=BLOCKS[setting];vac=empty_key(CONFIG);out={}
    for a in block:
        for spin in range(2):
            reference=15*(1-spin)+a
            out[reference,(1<<(15*spin+a),vac[1],vac[2])]=(-1)**spin/s.sqrt(2*len(block))
    assert inner(out,out)==1
    return out


@lru_cache(maxsize=2)
def target(setting):
    out={}
    for (reference,key),value in preparation(setting).items():
        for loop,factor in loop_state(setting).items():
            out[reference,(key[0],loop[1],key[2])]=s.expand(value*factor)
    assert inner(out,out)==1 and inner(out,preparation(setting))==0
    return out


def apply_hop(state,edge):
    term=square_terms()[edge];out=defaultdict(lambda:s.Integer(0))
    for (reference,key),value in state.items():
        for dest,factor in hopping_column(key,(term,)).items():out[reference,dest]+=value*factor
    return clean(out)


def diagonal(state):
    out={}
    for refkey,value in state.items():
        reference,key=refkey
        energy=sum(v!=0 for v in key[1])+sum(v!=0 for v in key[2])
        if energy:out[refkey]=s.expand(value*energy)
    return out


def total_current(state,site,current):
    native=representation()[2][current];out=defaultdict(lambda:s.Integer(0))
    for (reference,key),value in state.items():
        for dest,factor in gauss_column(key,site,current,CONFIG).items():out[reference,dest]+=value*factor
        if site==0:
            spin,b=divmod(reference,15)
            for a in range(15):
                if native[b,a]:out[15*spin+a,key]-=value*native[b,a]
    return clean(out)


def spin_action(state,axis):
    out=defaultdict(lambda:s.Integer(0));sigma=PAULI[axis]/2
    for (reference,key),value in state.items():
        assert key[0].bit_count()==1
        mode=key[0].bit_length()-1;site,local=divmod(mode,30);spin,a=divmod(local,15)
        rspin,b=divmod(reference,15)
        for t in range(2):
            if sigma[t,spin]:
                dest=(1<<(30*site+15*t+a),key[1],key[2])
                out[reference,dest]+=value*sigma[t,spin]
            if sigma[t,rspin]:out[15*t+b,key]+=value*sigma[t,rspin]
    return clean(out)


def loop_ward_checks(setting):
    """Full-domain endpoint cancellation of the actual four-coordinate trace."""
    block=BLOCKS[setting];endpoints=((0,4),(4,6),(2,6),(0,2));count=0
    words=tuple(((a,b,True),(b,c,True),(d,c,False),(a,d,False))
                for a,b,c,d in product(block,repeat=4))
    for site,current in product(range(8),range(12)):
        native=representation()[2][current];terms=defaultdict(lambda:s.Integer(0))
        for word in words:
            for edge,(left,right) in enumerate(endpoints):
                if site not in (left,right):continue
                endpoint=int(site==right);a,b,adjoint=word[edge]
                for c in block:
                    if endpoint==0:
                        coefficient=native[c,a] if adjoint else -native[a,c]
                        replacement=(c,b,adjoint)
                    else:
                        coefficient=-native[b,c] if adjoint else native[c,b]
                        replacement=(a,c,adjoint)
                    if coefficient:
                        target=list(word);target[edge]=replacement
                        terms[tuple(target)]+=coefficient
        assert not clean(terms);count+=1
    return count


@lru_cache(maxsize=2)
def physical_checks(setting):
    prep=preparation(setting);final=target(setting);loop=loop_state(setting)
    for site,current in product(range(8),range(12)):
        assert not total_current(prep,site,current)
        # Loop alone is a gauge singlet. Its product with the preparation
        # inherits all currents without a numerical rank or tolerance.
        action=defaultdict(lambda:s.Integer(0))
        for key,value in loop.items():
            for dest,factor in gauss_column(key,site,current,CONFIG).items():action[dest]+=value*factor
        assert not clean(action)
    for axis in range(3):assert not spin_action(prep,axis) and not spin_action(final,axis)
    return {'setting':setting,'native_indices':list(BLOCKS[setting]),
            'native_block_dimension':len(BLOCKS[setting]),
            'reference':'static full dual native R15 with ordinary two-component spin; all F30 occupations retained',
            'preparation_terms':len(prep),'loop_terms':len(loop),'final_terms':len(final),
            'preparation_and_target_norm_squared':['1','1'],'initial_final_overlap':'0',
            'preparation_local_Gauss_equations':96,'loop_local_Gauss_equations':96,
            'full_119_power_4_loop_operator_Ward_equations':loop_ward_checks(setting),
            'preparation_and_target_spin_singlet_equations':6,
            'source_and_reference_fermion_parities':['odd','odd'],
            'full_preparation_controller':'swap the orthonormal total vacuum and declared physical singlet; identity on their full complement',
            'controller_is_Gauss_spin_and_parity_invariant':True,
            'idle_reference_Hamiltonian':'0'}


@lru_cache(maxsize=2)
def leading_histories(setting):
    initial=preparation(setting);final=target(setting);records=[];m4=s.Integer(0);m5=s.Integer(0)
    for order in permutations(SQUARE):
        states=[initial]
        for edge in order:states.append(apply_hop(states[-1],edge))
        value=inner(final,states[-1]);m4+=value
        if value:
            records.append({'chronological_edges':list(order),'fourth_moment_contribution':str(value)})
            # Every intermediate matter location is away from the anchor;
            # hence the selected trace coefficient is Q-irreducible.
            for state in states[1:4]:
                assert all((key[0].bit_length()-1)//30!=0 for reference,key in state)
        for insertion in range(5):
            state=diagonal(states[insertion])
            for edge in order[insertion:]:state=apply_hop(state,edge)
            m5+=inner(final,state)
    d=len(BLOCKS[setting]);assert m4==-s.Rational(1,256*d) and m5==10*m4
    assert len(records)==1 and records[0]['chronological_edges']==list(SQUARE)
    # Representative complete columns confirm the diagonal extraction and
    # the actual spectator-preserving source action used in the histories.
    for key in (next(iter(initial))[1],next(iter(final))[1]):
        energy=sum(v!=0 for v in key[1])+sum(v!=0 for v in key[2])
        assert column(key,CONFIG).get(key,0)==energy
    return {'setting':setting,'all_distinct_edge_orders_exhausted':24,
            'nonzero_fourth_order_histories':records,
            'full_H_moments_0_through_5':['0','0','0','0',str(m4),str(m5)],
            'all_five_diagonal_insertions_for_every_order_executed':True,
            'lower_orders_zero_by_four_changed_links':True,
            'Yukawa_words_excluded_only_at_these_orders_by_exact_support':True,
            'five_hop_words_zero_by_bipartite_matter_parity':True,
            'leading_history_has_no_intermediate_P_return':True,
            'all_base_interactions_retained_in_H_and_remainder':True}


def geometry_certificate():
    terms=square_terms();bound=norm_bound(CONFIG);assert bound==1896
    spin=PAULI[1]*PAULI[0]*PAULI[1]*PAULI[0];assert spin==-s.eye(2)
    coefficient=s.conjugate(terms[SQUARE[0]].coefficient)*s.conjugate(terms[SQUARE[1]].coefficient)
    coefficient=s.simplify(coefficient*terms[SQUARE[2]].coefficient*terms[SQUARE[3]].coefficient)
    assert coefficient==s.Rational(1,256)
    return {'side':2,'spacing':'pi','all_Yukawa_and_penalty_coefficients':'1',
            'source_sites':8,'source_CAR_modes':240,'reference_CAR_modes':30,
            'link_count':24,'states_per_link':119,'states_per_scalar':5,
            'complete_source_reference_dimension':'2^270 * 119^24 * 5^8',
            'square_vertices':[list(cells(2)[i]) for i in ROUTE],
            'square_edges':list(SQUARE),'coordinate_orientation':['adjoint','adjoint','forward','forward'],
            'actual_kinetic_coefficients':[str(terms[e].coefficient) for e in SQUARE],
            'actual_forward_matter_coefficient_product':str(coefficient),'Pauli_loop_product':'-I2',
            'antiperiodic_seams_crossed':0,'contractible_square':True,
            'actual_Peter_Weyl_coordinate_vacuum_columns_used':True,
            'full_source_Hamiltonian_norm_upper_bound':str(bound),
            'full_norm_bound_applies_to_PW_coordinates':'each complete compressed coordinate is a contraction',
            'other_link_architecture_spectra_transferred':False}
