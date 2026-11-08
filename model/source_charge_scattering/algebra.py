"""Source charge blocks, pair preparation and exact matrix-algebra witnesses."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations,product
from math import comb
import numpy as np
from model.protected_defect_catalog import charge
from model.native_ordered_ground import source
from model.centered_particle_spectrum.certificate import decimal_bounds

PAIRS=tuple(combinations(range(8),2))


def local_creation(op,x,q):
    if type(x) is not int or not 0<=x<8 or q not in(-1,1):
        raise ValueError('actual cube site and signed unit charge required')
    identity=np.eye(256,dtype=np.int64);indices=np.arange(256)
    z=np.diag(1-2*((indices>>x)&1));xx=identity[indices^(1<<x)]
    return z@(identity-q*op['signs'][x]*xx)


def positive_leading_minors(matrix):
    """Fraction-free Sylvester certificate; all divisions are checked exactly."""
    a=[[int(x) for x in row] for row in matrix]
    n=len(a)
    if not n or any(len(row)!=n for row in a) or any(a[i][j]!=a[j][i] for i in range(n) for j in range(n)):
        raise ValueError('nonempty symmetric integer matrix required')
    previous=1;pivots=[]
    for k in range(n):
        pivot=a[k][k]
        if pivot<=0:raise ArithmeticError('positive Gram floor failed')
        pivots.append(pivot)
        if k==n-1:break
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=a[i][j]*pivot-a[i][k]*a[k][j]
                value,remainder=divmod(numerator,previous)
                if remainder:raise ArithmeticError('Bareiss exact division failed')
                a[i][j]=value
        for i in range(k+1,n):a[i][k]=0
        previous=pivot
    return tuple(pivots)


def spanning_walk(matrix):
    n=len(matrix);seen={0};walk=[0]
    def visit(i):
        for j in range(n):
            if j not in seen and matrix[i,j]!=0:
                seen.add(j);walk.append(j);visit(j);walk.append(i)
    visit(0)
    if len(seen)!=n:raise ArithmeticError('actual Hamiltonian graph is disconnected')
    if any(matrix[i,j]==0 for i,j in zip(walk,walk[1:])):
        raise ArithmeticError('walk uses a nonexistent source transition')
    return tuple(walk)


@lru_cache(maxsize=1)
def data():
    op=charge.operators();g,error,checked=charge.ground_input()
    g=np.asarray(g,dtype=object);norm=int(g@g);scale=op['scale']
    creation={(x,q):local_creation(op,x,q) for x in range(8) for q in(-1,1)}
    total=np.zeros((256,256),dtype=np.int64)
    for q,r in product((-1,1),repeat=2):
        k=creation[0,q]@creation[7,r]
        if not np.array_equal(op['twice_charge']@k-k@op['twice_charge'],2*(q+r)*k):
            raise ArithmeticError('two local creations must add their source charges')
        total+=k.T@k
    if not np.array_equal(total,16*np.eye(256,dtype=np.int64)):
        raise ArithmeticError('four complete separated creation outcomes required')
    sectors={}
    for q in(-1,1):
        labels=tuple(i for i,c in enumerate(op['charges']) if c==2*q)
        w=op['walsh'][:,labels].astype(object)
        if len(labels)!=28 or not np.array_equal(w.T@w,256*np.eye(28,dtype=object)):
            raise ArithmeticError('full charge-two basis must have28 orthogonal source labels')
        numerator=w.T@op['H'].astype(object)@w
        if any(int(v)%256 for v in numerator.flat):raise ArithmeticError('exact charge Hamiltonian denominator changed')
        block=numerator//256
        if not np.array_equal(op['H'].astype(object)@w,w@block):
            raise ArithmeticError('charge block is not a reducing part of the complete H')
        images=np.array([creation[x,q].astype(object)@creation[y,q].astype(object)@g for x,y in PAIRS],dtype=object).T
        coordinates=w.T@images
        den=4096*norm
        gram_num=coordinates.T@coordinates
        # Project only the untrusted approximation. True pair vectors have this
        # exact charge; all omitted proposal norm is bounded by the ground error.
        lost=[]
        for i in range(28):
            lost_i=F(int(images[:,i]@images[:,i]),16*norm)-F(int(gram_num[i,i]),den)
            if not 0<=lost_i<=error**2:raise ArithmeticError('proposal charge leakage exceeds source certificate')
            lost.append(lost_i)
        floor=F(1,10**9)
        shifted=gram_num*floor.denominator-den*floor.numerator*np.eye(28,dtype=object)
        pivots=positive_leading_minors(shifted)
        true_floor=floor-56*error
        if true_floor<=0:raise ArithmeticError('true source pair frame must be bounded below')
        walk=spanning_walk(block)
        sectors[q]={'labels':labels,'W':w,'block':block,'coordinates':coordinates,
            'denominator':den,'gram_numerator':gram_num,'pivots':pivots,
            'true_gram_floor':true_floor,'projected_proposal_leakage':tuple(lost),'walk':walk}
    return op,g,norm,error,checked,creation,sectors


def report():
    op,g,norm,error,checked,creation,sectors=data()
    records=[]
    for q,r in product((-1,1),repeat=2):
        v=creation[0,q].astype(object)@creation[7,r].astype(object)@g
        weight=F(int(v@v),16*norm)
        records.append({'first_charge':q,'second_charge':r,'total_charge':q+r,
            'original_weight_interval':decimal_bounds(max(F(),weight-2*error),min(F(1),weight+2*error),12)})
    return {'ground_source_certificate':checked,'pair_sites':[[x,y] for x,y in PAIRS],
        'complete_separated_input_pair':[0,7],'all_four_preparation_records':records,
        'same_charge_pair_sectors':[
            {'total_charge':2*q,'full_sector_dimension':28,'exact_X_labels':list(s['labels']),
             'source_H_integer_block':[[int(x) for x in row] for row in s['block']],
             'source_H_denominator':op['scale'],
             'source_nonzero_transition_walk':list(s['walk']),
             'all28_states_visited':len(set(s['walk']))==28,
             'proposal_Gram_floor':'1/1000000000',
             'true_selected_source_Gram_lower':str(s['true_gram_floor']),
             'Sylvester_integer_principal_minors':[str(x) for x in s['pivots']],
             'largest_proposal_weight_outside_exact_charge':str(max(s['projected_proposal_leakage'])),
             'all_pair_preparations_span_complete_reducing_sector':True,
             'diagonal_units_and_actual_H_generate_full_M28':True}
            for q,s in sectors.items()],
        'maximal_algebra_theorem':'OPH.QFT.adjoin_eq_top_of_walk',
        'application_kind':'analytic application to exact executable source inputs; no new Lean module compiled',
        'generic_native_instrument_compiler':'model/spin_operational_position/PROOF.md section3',
        'finite_dimensional_scaling_boundary':[
            {'sites':n,'two_local_pair_vectors':comb(n,2),'total_charge2_sector_dimension':comb(n,n//2+2)}
            for n in(8,10,12,16)],
        'full_two_particle_Fock_structure_assumed':False}
