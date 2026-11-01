"""Exact ordered-word reflection checks and Poisson tail bounds at fixed spacing.

The finite rational examples check ordering and signs. Existence of every
infinite-volume Euclidean n-point function follows from the positive GNS
generator and the dominated-convergence argument in README.md, not from
these examples. No rotational Euclidean or continuum axioms are inferred.
"""
from fractions import Fraction as F
from itertools import combinations


def _exact(value):
    if type(value) not in (int,F):
        raise ValueError('exact rational value required')
    return F(value)


def damping(energies,steps,base=F(1,2)):
    """exp(-steps*log(1/base)*L) for a nonnegative integer diagonal L."""
    if (not energies or any(type(e) is not int or e<0 for e in energies)
            or energies[0]!=0 or type(steps) is not int or steps<0):
        raise ValueError('nonnegative diagonal generator with stationary vacuum required')
    base=_exact(base)
    if not 0<base<=1:raise ValueError('contractive positive semigroup base required')
    return tuple(base**(steps*e) for e in energies)


def _apply(operator,vector):
    if len(operator)!=len(vector) or any(len(row)!=len(vector) for row in operator):
        raise ValueError('operator dimension changed')
    return tuple(sum((_exact(a)*b for a,b in zip(row,vector,strict=True)),F(0))
                 for row in operator)


def adjoint(operator):
    # These finite witnesses are real rational. The theorem uses the full
    # complex Hilbert-space adjoint, with no reality restriction.
    return tuple(tuple(row) for row in zip(*operator,strict=True))


def _word_domain(operators,times,energies,positive):
    damping(energies,0)
    if (not operators or len(operators)!=len(times)
            or any(type(t) is not int for t in times)
            or any(a>b for a,b in zip(times,times[1:]))
            or (positive and times[0]<0)):
        raise ValueError('nonempty ordered word and ordered integer times required')


def half_word(operators,times,energies):
    """e^(-t1 L) A1 ... e^(-(tn-tn-1)L) An Omega, in log(2) units."""
    _word_domain(operators,times,energies,True)
    vector=(F(1),)+(F(0),)*(len(energies)-1)
    for j in range(len(operators)-1,-1,-1):
        vector=_apply(operators[j],vector)
        gap=times[j]-(times[j-1] if j else 0)
        vector=tuple(a*b for a,b in zip(damping(energies,gap),vector,strict=True))
    return vector


def ordered_expectation(operators,times,energies):
    """Omega A1 e^(-(t2-t1)L) A2 ... An Omega; times may straddle zero."""
    _word_domain(operators,times,energies,False)
    vector=(F(1),)+(F(0),)*(len(energies)-1)
    for j in range(len(operators)-1,-1,-1):
        vector=_apply(operators[j],vector)
        if j:
            vector=tuple(a*b for a,b in zip(
                damping(energies,times[j]-times[j-1]),vector,strict=True))
    return vector[0]


def reflected_entry(left,right,energies):
    a,t=left;b,s=right
    _word_domain(a,t,energies,True);_word_domain(b,s,energies,True)
    operators=tuple(adjoint(x) for x in reversed(a))+tuple(b)
    times=tuple(-x for x in reversed(t))+tuple(s)
    return ordered_expectation(operators,times,energies)


def poisson_tail_bound(gaps,cutoffs):
    """Rational union bound using pi>3; zero gaps need no Poisson integral."""
    if len(gaps)!=len(cutoffs):raise ValueError('one cutoff per gap required')
    total=F(0)
    for gap,cutoff in zip(gaps,cutoffs,strict=True):
        gap,cutoff=_exact(gap),_exact(cutoff)
        if gap<0 or cutoff<=0:raise ValueError('nonnegative gap and positive cutoff required')
        total+=2*gap/(3*cutoff)
    return min(F(1),total)


def _annihilator(mode):
    rows=[[F(0) for _ in range(4)] for _ in range(4)]
    for mask in range(4):
        if mask>>mode&1:
            rows[mask^(1<<mode)][mask]=F((-1)**((mask&((1<<mode)-1)).bit_count()))
    return tuple(tuple(row) for row in rows)


def witness_words():
    a,b=_annihilator(0),_annihilator(1)
    aa,bb=adjoint(a),adjoint(b)
    return ((((aa,),(1,))),(((b,bb,aa),(1,2,3))),(((bb,),(1,))),
            (((aa,bb),(1,1))),(((bb,aa),(1,1))))


def _determinant(matrix):
    rows=[list(row) for row in matrix];answer=F(1)
    for j in range(len(rows)):
        pivot=next((i for i in range(j,len(rows)) if rows[i][j]),None)
        if pivot is None:return F(0)
        if pivot!=j:rows[pivot],rows[j]=rows[j],rows[pivot];answer=-answer
        value=rows[j][j];answer*=value
        for i in range(j+1,len(rows)):
            ratio=rows[i][j]/value
            rows[i]=[a-ratio*b for a,b in zip(rows[i],rows[j],strict=True)]
    return answer


def certificate():
    energies=(0,1,2,3);words=witness_words()
    vectors=[half_word(*word,energies) for word in words]
    expected=[(0,F(1,2),0,0),(0,F(1,32),0,0),(0,0,F(1,4),0),
              (0,0,0,F(1,8)),(0,0,0,-F(1,8))]
    if vectors!=expected:raise ArithmeticError('ordered CAR half-words changed')
    gram=tuple(tuple(reflected_entry(a,b,energies) for b in words) for a in words)
    if any(gram[i][j]!=sum(a*b for a,b in zip(vectors[i],vectors[j],strict=True))
           for i in range(5) for j in range(5)):
        raise ArithmeticError('reflection lost reversed adjoints or word order')
    minors=[_determinant([[gram[i][j] for j in indices] for i in indices])
            for size in range(1,6) for indices in combinations(range(5),size)]
    if min(minors)<0:raise ArithmeticError('reflected word Gram is not positive')
    coefficients=(1,-2,3,5,7)
    norm=sum(coefficients[i]*gram[i][j]*coefficients[j] for i in range(5) for j in range(5))
    if norm!=F(209,256) or gram[3][4]!=-F(1,64):
        raise ArithmeticError('fermionic reflection sign control changed')
    return {'all_ordered_local_Euclidean_npoint_limits_exist':True,
            'same_selected_ground_subsequence':True,
            'positive_time_reflection_is_word_Gram':True,
            'proof':'positive GNS generator; product Poisson representation; local real-time convergence and dominated convergence',
            'finite_witness_dimension':4,'reflected_words':5,
            'all_principal_minors_checked':len(minors),
            'exact_reflection_Gram':[[str(x) for x in row] for row in gram],
            'signed_polynomial_norm_squared':str(norm),
            'odd_word_cross_entry':str(gram[3][4]),
            'Poisson_tail_probability_bound_example':str(poisson_tail_bound((F(1,2),F(3,2)),(100,100))),
            'complex_time_boundary_rate_computed':False,
            'limiting_hierarchy_values_computed':False,
            'bosonic_permutation_invariance_inferred':False,
            'Euclidean_O4_covariance_inferred':False,
            'continuum_OS_axioms_completed':False}
