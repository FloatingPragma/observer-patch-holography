"""Finite attempt and local-repair-clock bounds from the literal source graph."""
from collections import defaultdict
from fractions import Fraction as F
from math import comb
from . import cone


def integer(value,minimum=0,maximum=None):
    if type(value) is not int or value<minimum or (maximum is not None and value>maximum):
        raise ValueError('an exact integer in the stated finite evaluator domain required')
    return value


def ball(family,root,radius):
    integer(radius);cone.subset(root,len(family.weights));seen=set(root)
    for _ in range(radius):
        seen |= {v for edge in family.seams if seen.intersection(edge) for v in edge}
    return tuple(sorted(seen))


def escape(family,root,window,attempts,*,max_states=65536):
    """Exact selector-only cone probability with an absorbing escape flag."""
    integer(attempts,0,256);integer(max_states,1)
    cone.subset(root,len(family.weights));cone.subset(window,len(family.weights))
    if not root or not set(root)<=set(window):
        raise ValueError('nonempty roots must lie in the comparison window')
    states={frozenset(root):F(1)};bad=F();largest=1
    selected=tuple((edge,F(weight,sum(family.seam_weights)))
                   for edge,weight in zip(family.seams,family.seam_weights,strict=True))
    for _ in range(attempts):
        following=defaultdict(F)
        for support,probability in states.items():
            for edge,selection in selected:
                target=support|set(edge) if support.intersection(edge) else support
                mass=probability*selection
                if target<=set(window):following[frozenset(target)]+=mass
                else:bad+=mass
        if len(following)>max_states:
            raise ValueError('exact cone subset evaluator resource stop')
        states=dict(following);largest=max(largest,len(states))
    if bad+sum(states.values())!=1:
        raise ArithmeticError('cone selector probabilities do not normalize')
    return {'escape_probability':bad,'surviving_subsets':states,'maximum_subset_states':largest}


def walk_bound(family,root,window,attempts):
    """Union over time-ordered first-exit walks, including every walk length.

    The finite graph vector recursion avoids full endpoint-state matrices.
    It is a conservative walk bound; the exact cone can be smaller.
    """
    integer(attempts,0,64);cone.subset(root,len(family.weights));cone.subset(window,len(family.weights))
    if not root or not set(root)<=set(window):raise ValueError('roots must lie in the window')
    inside=set(window);vector={v:F(int(v in root)) for v in window};coefficients=[]
    selected=tuple((edge,F(weight,sum(family.seam_weights)))
                   for edge,weight in zip(family.seams,family.seam_weights,strict=True))
    for _ in range(attempts):
        following=defaultdict(F);coefficient=F()
        for (a,b),q in selected:
            if a in inside and b in inside:
                following[b]+=vector[a]*q;following[a]+=vector[b]*q
            elif a in inside:coefficient+=vector[a]*q
            elif b in inside:coefficient+=vector[b]*q
        coefficients.append(coefficient)
        vector={v:following[v] for v in window}
    total=sum((comb(attempts,k)*x for k,x in enumerate(coefficients,1)),F())
    return {'first_exit_walk_coefficients':tuple(coefficients),'uncapped_union_bound':total,
            'error_upper':min(F(1),total)}


def radius_bound(family,root,radius,attempts):
    integer(radius);integer(attempts)
    cone.subset(root,len(family.weights))
    degrees=defaultdict(int)
    for a,b in family.seams:degrees[a]+=1;degrees[b]+=1
    length=radius+1
    count=comb(attempts,length) if length<=attempts else 0
    q=F(max(family.seam_weights),sum(family.seam_weights))
    return min(F(1),len(root)*count*(max(degrees.values())*q)**length)


def halo_geometry(level,cell):
    graph=cone.native.geometry.federation(level);integer(cell,0,graph.cells-1)
    root=tuple(range(12*cell,12*(cell+1)));root_set=set(root)
    halo=tuple(sorted({v for u in root for v in graph.neighbours[u]}-root_set))
    window=root_set|set(halo)
    if (len(halo)!=3 or len({v//12 for v in halo})!=3
            or any(set(graph.neighbours[v])&set(halo) for v in halo)
            or any(len(set(graph.neighbours[v])&root_set)!=1 or
                   len(set(graph.neighbours[v])-window)!=5 for v in halo)):
        raise ArithmeticError('the literal three-port halo does not have the required ownership')
    if sum(set(e)<=root_set for e in graph.seams)!=30:
        raise ArithmeticError('the native cell lost its thirty internal seams')
    return {'root':root,'halo':halo,'window':tuple(sorted(window)),
            'seams':len(graph.seams),'internal_seams':30,'birth_seams':3,
            'outside_seams_per_halo_port':5,'ordered_two_edge_patterns':15}


def matrix_apply(matrix,vector):
    return tuple(sum((a*b for a,b in zip(row,vector,strict=True)),F()) for row in matrix)


def halo_transition(seams):
    integer(seams,45)
    matrix=[[F() for _ in range(4)] for _ in range(4)]
    for n in range(4):
        matrix[n][n]=1-F(3+4*n,seams)
        if n<3:matrix[n+1][n]=F(3-n,seams)
    return tuple(map(tuple,matrix))


def halo_global_error(seams,attempts):
    integer(attempts,0,256);matrix=halo_transition(seams);vector=(F(1),F(),F(),F())
    for _ in range(attempts):vector=matrix_apply(matrix,vector)
    error=1-sum(vector)
    upper=min(F(1),F(15*attempts*(attempts-1),2*seams*seams))
    if not 0<=error<=upper:
        raise ArithmeticError('exact halo escape exceeds the ordered-pair union bound')
    return {'cone_escape_probability':error,'complete_local_halfdiamond_upper':error,
            'ordered_pair_upper':upper,'no_raw_output_distance_equality_asserted':True}


def internal_block():
    """30(D-B)^-1, the exact reversed geometric failure-block kernel."""
    result=[[F() for _ in range(4)] for _ in range(4)]
    for start in range(4):
        value=F(30,33+4*start);result[start][start]=value
        for end in range(start+1,4):
            value*=F(4-end,33+4*end);result[end][start]=value
    return tuple(map(tuple,result))


def halo_internal_error(repairs):
    integer(repairs,0,256);matrix=internal_block();vector=(F(1),F(),F(),F())
    for _ in range(repairs):vector=matrix_apply(matrix,vector)
    error=1-sum(vector);upper=min(F(1),F(repairs*(repairs+1),120))
    if not 0<=error<=upper:
        raise ArithmeticError('geometric-block escape exceeds the ordered-pair stopping bound')
    return {'internal_repairs':repairs,'cone_escape_probability':error,
            'complete_local_halfdiamond_upper':error,'ordered_pair_upper':upper,
            'global_elapsed_attempt_record_retained':True,
            'cone_escape_bound_independent_of_global_seam_count':True,
            'stops_immediately_at_final_internal_repair':True}


def fixed_timeout_bound(family,root,window,repairs,attempts):
    """Generic local-clock comparison with an exact finite binomial timeout."""
    integer(repairs,1);integer(attempts,0,256)
    inside=set(root)
    rate=sum((family.selection(e) for e in family.seams if set(e)<=inside),F())
    if not rate:raise ValueError('the chosen local repair clock has no internal seam')
    timeout=sum((comb(attempts,j)*rate**j*(1-rate)**(attempts-j)
                 for j in range(min(repairs,attempts+1))),F())
    cone_error=escape(family,root,window,attempts)['escape_probability']
    return {'internal_attempt_probability':rate,'timeout_probability':timeout,
            'cone_error_at_timeout':cone_error,'error_upper':min(F(1),timeout+cone_error)}


def abel_response_bound(q):
    """Exact geometric averaging of bounded local multi-time record scores."""
    if type(q) is not F or not 0<=q<1:
        raise ValueError('a declared rational analysis parameter0<=q<1 required')
    matrix=internal_block();solution=[]
    for i in range(4):
        rhs=F(i==0)+q*sum((matrix[i][j]*solution[j] for j in range(i)),F())
        solution.append(rhs/(1-q*matrix[i][i]))
    bound=2*(1-(1-q)*sum(solution))
    coarse=min(F(2),q/(30*(1-q)**2))
    if not 0<=bound<=coarse:
        raise ArithmeticError('exact Abel-response error exceeds the ordered-pair envelope')
    return {'analysis_parameter':q,'bounded_score_response_error_upper':bound,
            'coarse_error_upper':coarse,'resolvent_vector':tuple(solution),
            'infinite_series_is_an_analysis_of_the_kernel_family':True,
            'unbounded_native_program_execution_claimed':False,
            'physical_energy_or_scattering_parameter':False}
