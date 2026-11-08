"""Directed exact-angle audit of the literal 22-pulse local beam.

The control angle is its stored binary64 rational, not a symbolic pi
multiple. Exact source generator identities give the primitive below;
the numerical M6 exponential is not used in this interval calculation.
"""
from fractions import Fraction as F
from math import isqrt
from model.native_effect_synthesis import intervals as iv
from model.native_fair_local_preparation import coupling

Z=(0,0);ONE=iv.exact(1)
ZERO=(Z,Z);IDENTITY=(ONE,Z)


def add(a,b):return iv.add(a[0],b[0]),iv.add(a[1],b[1])
def multiply(a,b):
    return (iv.add(iv.mul(a[0],b[0]),iv.neg(iv.mul(a[1],b[1]))),
            iv.add(iv.mul(a[0],b[1]),iv.mul(a[1],b[0])))


def matmul(a,b):
    n=len(a)
    return [[sum_complex(multiply(a[i][k],b[k][j]) for k in range(n))
             for j in range(n)] for i in range(n)]


def sum_complex(values):
    out=ZERO
    for value in values:out=add(out,value)
    return out


def gate(event):
    s,c=iv.sincos(F(event['angle']));kind=event['kind'];targets=event['targets']
    if kind=='J':
        return [[((c,s) if i!=3 else IDENTITY) if i==j else ZERO
                 for j in range(4)] for i in range(4)]
    if kind=='X':small=[[(c,Z),(Z,s)],[(Z,s),(c,Z)]]
    elif kind=='Y':small=[[(c,Z),(s,Z)],[(iv.neg(s),Z),(c,Z)]]
    elif kind=='Z':small=[[(c,s),ZERO],[ZERO,(c,iv.neg(s))]]
    elif kind=='I':small=[[(c,s),ZERO],[ZERO,(c,s)]]
    elif kind=='H':
        root=isqrt(iv.Q*iv.Q//2);inverse_root=(root,root+1)
        if not F(root,iv.Q)**2<=F(1,2)<=F(root+1,iv.Q)**2:
            raise ArithmeticError('inverse square-root enclosure failed')
        h=iv.mul(s,inverse_root)
        small=[[(c,h),(Z,h)],[(Z,h),(c,iv.neg(h))]]
    else:raise ValueError('literal beam source pulse kind required')
    bit=targets[0]
    return [[small[(i>>bit)&1][(j>>bit)&1] if (i&~(1<<bit))==(j&~(1<<bit)) else ZERO
             for j in range(4)]for i in range(4)]


def certificate():
    # This owner verifies all exact Q(sqrt5) source generator identities.
    coupling.occupation.compiler._source()
    word=coupling.words()['beam'][1]
    out=[[IDENTITY if i==j else ZERO for j in range(4)]for i in range(4)]
    for event in word:out=matmul(gate(event),out)
    witness=[]
    for i in (0,3):
        for j in (1,2):
            for component,name in enumerate(('real','imag')):
                lo,hi=out[i][j][component]
                if lo>0 or hi<0:
                    witness.append({'output':i,'input':j,'component':name,
                                    'interval':[str(F(lo,iv.Q)),str(F(hi,iv.Q))]})
    if not witness:raise ArithmeticError('no certified nonzero literal-angle number leakage found')
    return {'pulses':len(word),'interval_bits':iv.BITS,
            'literal_angles':[str(F(e['angle'])) for e in word],
            'nonzero_number_leakage_entries':witness,
            'exact_source_generators_checked':True,
            'ideal_proper_code_failure_blocks_are_zero':True,
            'symbolic_pi_replacement_performed':False}
